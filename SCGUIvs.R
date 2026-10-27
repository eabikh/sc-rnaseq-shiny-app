# =============================================================================
# Single-Cell RNA-seq Analysis App  —  Seurat v5 | 1 / 2 / N samples
#
# QC, integration, and clustering, followed by cell-type annotation (ScType
# automated + manual override) and comparative analysis between conditions,
# with pathway-level interpretation of the resulting differential expression.
#
# Statistical methods:
#  - DE: Seurat::FindMarkers() (test.use as currently configured in the UI)
#  - GSEA: fgseaMultilevel() over msigdbr gene sets, ranked by
#    sign(logFC)*-log10(p)*abs(logFC)
#  - ORA: clusterProfiler::enricher() on the significant up/down genes from
#    a DE comparison
#
# Major feature areas:
#  - QC & clustering (filtering, PCA, integration, dimensionality reduction)
#  - ScType + manual cluster annotation
#  - Comparative DE with cached results (rv$de_cache), supporting multiple
#    simultaneous comparisons without re-running
#  - GSEA + Pathway Explorer (per-pathway bar/trace views)
#  - ORA (combined up/down plot, colored by direction)
#  - Global vs. per-cluster scope toggle for comparative violin/dot plots
#
# All plot downloads render at 600 DPI with dynamic, content-driven sizing
# (no fixed/clipped dimensions).
# =============================================================================

suppressPackageStartupMessages({
  library(shiny); library(Seurat); library(ggplot2); library(dplyr)
  library(patchwork); library(DT); library(writexl); library(hdf5r)
  library(ggrepel); library(fgsea); library(msigdbr); library(clusterProfiler)
})
options(shiny.maxRequestSize = 500 * 1024^2)

# =============================================================================
# HELPERS
# =============================================================================
`%||%` <- function(a,b) if(!is.null(a)) a else b

nat_sort <- function(x) {
  has_prefix <- grepl("^[0-9]+", x)
  if(all(has_prefix)) {
    num_prefix <- as.numeric(sub("^([0-9]+).*", "\\1", x))
    suffix <- sub("^[0-9]+", "", x)
    x[order(num_prefix, suffix)]
  } else sort(x)
}

# Fixed, ordered, stable color scale for identity/sample coloring — values
# must always come from the full, unsubsetted object so colors/order stay
# identical regardless of which cells are actually shown in a given view.
fixed_ident_scale <- function(values){
  lv <- nat_sort(unique(as.character(values)))
  ggplot2::scale_color_manual(values=setNames(scales::hue_pal()(length(lv)), lv), breaks=lv, drop=FALSE)
}

# Same color-per-level mapping as fixed_ident_scale, but for the `fill`
# aesthetic — UMAP points use `color`, violin shapes use `fill`.
fixed_ident_fill_scale <- function(values){
  lv <- nat_sort(unique(as.character(values)))
  ggplot2::scale_fill_manual(values=setNames(scales::hue_pal()(length(lv)), lv), breaks=lv, drop=FALSE)
}

get_mt_pat <- function(sp) switch(sp,"human"="^MT-","mouse"="^mt-","^MT-|^mt-|^Mt-")


apply_ident <- function(obj, ch) {
  if(ch=="sctype"&&"sctype_classification"%in%colnames(obj@meta.data))
    Idents(obj)<-obj$sctype_classification
  else if(ch=="manual"&&"manual_annotation"%in%colnames(obj@meta.data))
    Idents(obj)<-obj$manual_annotation
  else if(ch=="subcluster"&&"subcluster_labels"%in%colnames(obj@meta.data))
    Idents(obj)<-obj$subcluster_labels
  else Idents(obj)<-obj$seurat_clusters
  obj
}

resolve_genes <- function(raw, obj) {
  ag<-rownames(obj); found<-c()
  for(g in trimws(unlist(strsplit(raw,",")))) {
    if(!nchar(g)) next
    m<-grep(paste0("^",g,"$"),ag,ignore.case=TRUE,value=TRUE)
    if(!length(m)) m<-grep(g,ag,ignore.case=TRUE,value=TRUE)
    if(length(m)) found<-c(found,m[1])
  }; unique(found)
}

get_id_col <- function(ch) switch(ch%||%"cluster","sctype"="sctype_classification","manual"="manual_annotation","subcluster"="subcluster_labels","seurat_clusters")

# AggregateExpression() sanitises its group.by values inside agg@meta.data
# itself (verified empirically) — e.g. a purely-numeric identity like "3"
# is stored back as "g3", and "1_a" as "g1-a" — so the resulting column
# can't be found by string-matching the raw identity/sample values against
# agg@meta.data, and can't be reconstructed by guessing the sanitisation
# rule either. What IS preserved is the RELATIVE ORDER AggregateExpression
# groups (identity, sample) combinations in — nested by the identity
# column's own factor level order (or alphabetical if plain character),
# then by sample — with missing combinations simply omitted in place. So
# rank the raw (id, sample) combination against every other existing
# combination in that same order, computed straight from the original
# object, and use that rank to index into agg's own column names.
find_pb_col <- function(agg, id_col_val, sample_val, id_col_name, orig_obj){
  id_series <- orig_obj@meta.data[[id_col_name]]
  sample_series <- orig_obj@meta.data$sample
  combo_tab <- table(id_series, sample_series)
  present <- which(combo_tab > 0, arr.ind = TRUE)
  present <- present[order(present[,1], present[,2]), , drop=FALSE]
  ids   <- rownames(combo_tab)[present[,1]]
  samps <- colnames(combo_tab)[present[,2]]
  idx <- which(ids==as.character(id_col_val) & samps==as.character(sample_val))
  if(length(idx)) colnames(agg)[idx[1]] else NA_character_
}

# Sidebar "Active Samples" filter — display-only, never touches
# statistics/exports/composition. No-op if selection is empty or matches
# every currently-loaded sample, so an accidental "uncheck everything"
# never produces a blank plot.
filter_by_samples <- function(obj, active){
  if(is.null(active) || !length(active)) return(obj)
  cur <- unique(as.character(obj$sample))
  if(setequal(active, cur)) return(obj)
  subset(obj, subset = sample %in% active)
}

dl600 <- function(pfn,fname,w=12,h=8)
  downloadHandler(filename=function() paste0(fname,"_",Sys.Date(),".png"),
                  content=function(f){p<-pfn();ggsave(f,p,width=w,height=h,dpi=600,bg="white")})

# ---- cosmetic layout helpers (sizing only, no statistical effect) ----
px_from_in <- function(inches,ppi=92) round(inches*ppi)
max_label_len <- function(labels) if(length(labels)) max(nchar(as.character(labels)),0) else 0

# Split DotPlot (split.by) — width grows with #features, height with #rows (id x split).
# The 10in/8in floors below are only appropriate once there's enough
# features/rows to actually fill that canvas; below ~4 they're scaled down
# so a near-empty plot (e.g. 1 feature x 1 identity) gets a canvas sized to
# match, instead of leaving a handful of dots stranded in a mostly-blank
# 10x8in image. A sqrt ramp (not linear) is used so the floor doesn't
# shrink all the way toward 0 -- a straight linear ramp was tried first and
# measured to clip the split-violin's y-axis identity labels and its fixed
# with_dot_legend() legend column, both of which don't shrink with feature/
# row count; sqrt keeps a more generous floor at n=1 (50% of the full floor
# vs linear's 25%) while still tapering to the same 10in/8in at n>=4.
# dot.scale itself is never touched here -- this only resizes the canvas.
dotplot_dims <- function(n_features,n_rows,labels=character(0)){
  ml<-max_label_len(labels)
  w_floor<-10*sqrt(pmin(1,n_features/4))
  h_floor<-8*sqrt(pmin(1,n_rows/4))
  list(w=max(w_floor,n_features*1.2,ml*0.08), h=max(h_floor,n_rows*0.45))
}

# Same idea as cg_dot_scale, tuned to the Gene Expression tab's own dims
# formula (identities only, no split.by) and its own default dot.scale.
gene_dot_scale <- function(dims,n_features,n_rows,base=6){
  base
}
# coord_flip composition bar — height grows with #categories, width grows
# with label length AND keeps a real minimum panel width after the label
# margin is subtracted, so bars/dots always have genuine room to spread out.
# Previously this floored at a flat 10in regardless of label length; with
# long pathway descriptions, the matching plot.margin + auto-sized axis-text
# area (see ora_make_plot/gsea_make_plot's l=max(80,max_label_len*5) margin)
# can eat most of that 10in canvas, leaving only a sliver for the actual
# panel — measured empirically: a 55-char label left ~1.9in of usable panel
# out of a 10in canvas. The formula below (~8in baseline + 0.15in/label-char)
# was fit against measured panel width so a ~6in panel plus legend headroom
# survives regardless of label length; for short labels it still floors at 10.
barflip_dims <- function(n_cat,labels=character(0)){
  ml<-max_label_len(labels)
  list(w=max(10, 8+0.15*ml), h=max(8,n_cat*0.5))
}
# dodged composition bar (rotated x labels) — width grows with #categories
bardodge_dims <- function(n_cat,labels=character(0)){
  list(w=max(10,n_cat*0.8), h=max(7,max_label_len(labels)*0.06))
}
# DimPlot with repel labels — width grows modestly with label length so
# ggrepel has room; height is independent (forcing w=h distorted the point
# cloud since umap_1/umap_2 coordinate ranges aren't actually square).
repel_dims <- function(labels,base_w=12,base_h=8){
  ml<-max_label_len(labels)
  list(w=max(base_w,base_w+max(0,ml-15)*0.12), h=base_h)
}
repel_margin <- theme(plot.margin=margin(t=15,r=40,b=15,l=15,unit="pt"))

# Fixed minimum/maximum font sizes for legend + axis-title text so they stay
# legible regardless of how large the surrounding panel grows (violin/dot/UMAP).
min_text_theme <- theme(
  legend.text=element_text(size=9), legend.title=element_text(size=10),
  axis.title=element_text(size=10), axis.text=element_text(size=9),
  legend.key.height = unit(0.42,"cm"), legend.key.width = unit(0.42,"cm")
)

# Single source of truth for UMAP plot construction — used by every UMAP
# renderPlot and downloadHandler so on-screen and downloaded versions of
# the same plot never drift apart (theme/points/labels always match).
umap_plot <- function(obj, group.by = NULL, split.by = NULL, title = NULL,
                       label = FALSE, pt.size = 0.4, ncol = NULL,
                       reduction = "umap") {
  p <- Seurat::DimPlot(obj, reduction = reduction, group.by = group.by,
                        split.by = split.by, label = label, repel = label,
                        pt.size = pt.size, ncol = ncol)
  if (!is.null(title)) p <- p + ggplot2::ggtitle(title)
  p + repel_margin + ggplot2::coord_cartesian(clip = "off") + min_text_theme
}

# Full UMAP embedding's x/y range for the given object — used to pin
# FeaturePlot's axis range to match every other UMAP in the app, without
# touching panel/container size (coord_cartesian is a pure axis-range
# overlay, safe to apply even to a multi-gene FeaturePlot patchwork).
umap_full_xylim <- function(obj){
  emb <- Embeddings(obj, "umap")
  list(x=range(emb[,1]), y=range(emb[,2]))
}

# Fixed absolute panel size for UMAP plots, so the data panel occupies the
# same visual footprint regardless of label/legend width — renderPlot()
# does NOT auto-draw the resulting gtable (only ggplot objects), so every
# renderPlot call site must end with grid::grid.draw(fix_panel_size(...));
# ggsave()/dl600() already handle gtable objects directly, no grid.draw needed.
UMAP_PANEL_W_IN <- 8.5
UMAP_PANEL_H_IN <- 6.5
fix_panel_size <- function(p, w_in = UMAP_PANEL_W_IN, h_in = UMAP_PANEL_H_IN) {
  g   <- ggplot2::ggplotGrob(p)
  pan <- g$layout[g$layout$name == "panel", ]
  g$widths[pan$l]  <- grid::unit(w_in, "in")
  g$heights[pan$t] <- grid::unit(h_in, "in")
  g
}

# Inverse-scale a font size with item count: shrinks as n grows, never below
# `floor_v`, never above `ceiling_v` (used for marker-heatmap row/col labels).
mkr_font <- function(n,floor_v,ceiling_v,k) max(floor_v,min(ceiling_v,k/max(n,1)))

# Marker heatmap: gene/cluster label size scales inversely with count
# (floor so it never goes illegible, ceiling so it doesn't balloon with few
# clusters); top cluster color-key strip pinned to a small fixed proportion
# of plot height instead of ballooning with label size.
build_mkr_heatmap <- function(obj,genes,nc,ng,group_col="seurat_clusters"){
  y_sz<-mkr_font(ng,6,9,150); x_sz<-mkr_font(nc,7,10,80)
  obj@meta.data[[group_col]] <- factor(as.character(obj@meta.data[[group_col]]), levels=nat_sort(unique(as.character(obj@meta.data[[group_col]]))))
  DoHeatmap(obj,features=genes,group.by=group_col,group.bar.height=0.01,size=3)+NoLegend()+
    theme(axis.text.y=element_text(size=y_sz),axis.text.x=element_text(size=x_sz),
          plot.margin=margin(t=10,r=max(60, max_label_len(levels(obj@meta.data[[group_col]]))*4),b=10,l=10,unit="pt"))
}

# Comparative Gene Expression global/per-cluster toggle — when scope is
# "percluster" and a cluster is chosen, subsets to that identity only;
# otherwise returns obj unchanged (global behavior, today's default).
cg_scope_obj <- function(obj,id_col,scope,clust){
  if(identical(scope,"percluster") && !is.null(clust) && nchar(clust)){
    Idents(obj) <- obj@meta.data[[id_col]]
    obj <- subset(obj,idents=clust)
  } else if(identical(scope,"selected") && !is.null(clust) && length(clust)){
    Idents(obj) <- obj@meta.data[[id_col]]
    obj <- subset(obj,idents=clust)
  }
  obj
}

# Seurat's VlnPlot(stack=TRUE,flip=TRUE) internally hits rep(x,times=length(x)-2)
# when building its abbreviated y-axis tick labels; if a feature is perfectly
# constant within the plotted cells (common after subsetting to one small
# cluster, where a marker gene from another cluster reads as flat zero), the
# axis-break count collapses below 2 and that rep() gets a negative `times`,
# throwing "invalid 'times' argument". Drop constant features before plotting.
nonconstant_genes <- function(obj,genes){
  present<-genes[genes%in%rownames(obj)]
  if(!length(present)) return(present)
  d<-GetAssayData(obj,layer="data")[present,,drop=FALSE]
  keep<-apply(d,1,function(x) length(unique(x))>1)
  present[keep]
}

# Bottom-margin (and matching extra container height) for the rotated
# identity tick labels on cg_vln's x-axis (VlnPlot flip=TRUE puts identity
# there, not on the left) — sized off max(nchar(labels)), never off how many
# labels there are, so a single long per-cluster label still gets full room.
cg_vln_bottom_margin <- function(id_labels)
  max(90, max_label_len(id_labels) * 6)

# Extra left/right canvas room for long rotated identity labels.
# The bottom margin gives labels vertical room, but rotated labels can also
# extend beyond the left or right edge of the plot.
cg_vln_side_margin <- function(id_labels)
  max(40, max_label_len(id_labels) * 4)

# Shared width calculation for Comparative violin plots.
# Unlike the old formula, this does not stop growing at 14 inches when many
# identities or particularly long identity names are shown.
#
# cg_vln_min_canvas_in below adds a floor that guarantees a real panel
# survives once cg_vln_side_margin/cg_vln_bottom_margin's rotated-label
# margins are carved out of the canvas. Those margins (and ggplot's own
# auto-sized rotated-tick-label grob, which consumes *additional* space
# beyond the explicit margin) grow with label length, while the terms
# below don't grow at all for a single (or few) categories -- so a long
# identity label could consume nearly the entire canvas in both
# dimensions and leave a near-invisible sliver panel. Measured directly
# (render to PNG, detect the panel's actual pixel bounding box) with a
# single-category split-violin plot at several label lengths:
#   nchar=2  -> panel needs canvas w>=5.5in,  h>=7.06in to stay visible
#   nchar=32 -> panel needs canvas w>=9.87in, h>=11.53in
# The linear fits below (slope from those two measured points) reproduce
# that and keep a >=2in panel in both dimensions at any label length in
# between (verified at nchar=18 too); they floor at today's values for
# nchar<=2 so short-label plots are pixel-identical to before.
cg_vln_min_canvas_in <- function(id_labels) {
  ml <- max_label_len(id_labels)
  list(
    w = 5.5 + 0.162 * pmax(0, ml - 2),
    h = 7.06 + 0.1487 * pmax(0, ml - 2)
  )
}

cg_vln_width_in <- function(id_labels)
  max(
    5.5,
    length(id_labels) * 0.9,
    max_label_len(id_labels) * 0.15,
    cg_vln_min_canvas_in(id_labels)$w
  )

# Manual color-scale legend for split DotPlot (Seurat's DotPlot drops the
# Average Expression colorbar whenever split.by is combined with multiple
# plain `cols`, since it switches to scale_color_identity() internally).
# Builds one grey->sample-color gradient swatch per sample, each a compact
# vertical colorbar (gradient runs top-to-bottom, standard guide_colorbar
# convention) pinned to a small fixed size — not a stretched horizontal swatch.
dotplot_legend_bar <- function(sample_name,col){
  df<-data.frame(y=seq(0,1,length.out=50),x=1)
  bar<-ggplot(df,aes(x,y,fill=y))+geom_tile()+
    scale_fill_gradient(low="grey90",high=col,guide="none")+
    scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
    theme_void()+theme(plot.margin=margin(0,0,0,0,unit="pt"))
  title_grob<-ggplot()+theme_void()+labs(title=sample_name)+
    theme(plot.title=element_text(size=9,hjust=0.5,face="bold"),
          plot.margin=margin(2,2,2,2,unit="pt"))
  narrow_bar<-wrap_plots(list(plot_spacer(),bar,plot_spacer()),nrow=1,
                         widths=c(1,0.6,1))
  wrap_plots(list(title_grob,narrow_bar),ncol=1,heights=c(0.5,1))
}
dotplot_color_legend <- function(sample_names,cols,bar_cm=2.2){
  n<-length(sample_names); cols<-cols[seq_len(n)]
  bars<-lapply(seq_len(n),function(i) dotplot_legend_bar(sample_names[i],cols[i]))
  # fixed-cm height per swatch + a "null" spacer absorbing leftover column
  # space, so the colorbars stay compact at top instead of stretching to
  # fill however tall the main panel happens to be
  wrap_plots(c(bars,list(plot_spacer())),ncol=1)+
    plot_layout(heights=grid::unit(c(rep(bar_cm,n),1),c(rep("cm",n),"null")))+
    plot_annotation(title="Avg Expression\n(low → high)",
                     theme=theme(plot.title=element_text(size=9,hjust=0.5)))
}
with_dot_legend <- function(main_plot,sample_names,cols,legend_cm=1.8)
  (main_plot + dotplot_color_legend(sample_names,cols) +
     plot_layout(widths=grid::unit(c(1,legend_cm),c("null","cm"))))

# Per-cluster (single-category) styling for cg_vln only — Global mode is
# untouched since this is only invoked when scope=="percluster":
#  - Issue 1: Seurat's VlnPlot hardcodes scale="width"/trim=TRUE with a
#    zero-expansion y-scale internally; with only one category the kernel
#    density estimate fills out to the data edges much more (fewer points,
#    less taper), so peaked genes visually touch/clip the facet border.
#    `scale=`/`trim=` aren't exposed via VlnPlot's args, so the fix is to
#    re-add y-axis expansion afterward (verified empirically to resolve it).
#  - Issue 3: with only one x-axis category there's no overlap risk, so the
#    identity label can render larger and horizontal instead of rotated/small.
#  - Issue A: geom_violin(scale="width") normalizes the violin's width to
#    fill its *category slot*, and with only 1 real category, that one slot
#    is given the entire panel's discrete x-range, expanding it dramatically.
#    `scale=` isn't exposed by VlnPlot, but the slot width is purely a
#    function of how many positions exist on the discrete x scale — so we
#    add invisible/unlabeled "padding" categories around the real one. This
#    doesn't touch VlnPlot's internals or the data, only the already-built
#    ggplot's x scale, and only when is_single (never affects Global mode).
# Fixed canvas width (inches) for cg_vln when is_single — the padding trick
# above fixes the violin's shape/proportion *within* the panel, but the
# overall plot width was still sized by the multi-cluster formula, leaving
# dead space around the now-narrow violin. This shrinks the canvas itself.
cg_vln_single_width_in <- 5.5

cg_vln_apply_single <- function(p,is_single,real_id){
  if(!is_single) return(p)
  pads<-paste0("__cgvln_pad",1:4)
  lims<-c(pads[1],pads[2],real_id,pads[3],pads[4])
  p + scale_x_discrete(limits=lims,breaks=real_id) +
    scale_y_continuous(expand=expansion(mult=c(0.05,0.12))) +
    theme(axis.text.x=element_text(angle=0,hjust=0.5,size=12))
}

# =============================================================================
# GSEA helpers — reuse existing sizing/legend/font conventions, don't reinvent
# =============================================================================

# Ranking formula ported from bbulk.R exactly: sign(logFC)*-log10(pvals)*abs(logFC),
# with the same min_nonzero p-value floor for p=0 cases. Column names adapted
# for Seurat::FindMarkers output (gene/p_val/avg_log2FC instead of Gene/PValue/logFC).
build_gsea_stats <- function(de_res){
  pvals<-de_res$p_val
  if(all(is.na(pvals))) stats<-rep(NA_real_,nrow(de_res))
  else{
    min_nonzero<-min(pvals[pvals>0&is.finite(pvals)],na.rm=TRUE)
    pvals[!is.finite(pvals)|pvals<=0]<-min_nonzero
    stats<-sign(de_res$avg_log2FC)*-log10(pvals)*abs(de_res$avg_log2FC)
  }
  names(stats)<-de_res$gene
  stats<-stats[is.finite(stats)]
  sort(stats,decreasing=TRUE)
}

# Single continuous-gradient vertical colorbar (e.g. -log10(padj)), matching
# with_dot_legend()'s architecture: pinned-cm width column, compact fixed-cm
# swatch + null spacer so it stays a small key, not a stretched legend.
gsea_legend_bar <- function(low_col,high_col,title){
  df<-data.frame(y=seq(0,1,length.out=50),x=1)
  bar<-ggplot(df,aes(x,y,fill=y))+geom_tile()+
    scale_fill_gradient(low=low_col,high=high_col,guide="none")+
    scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+
    theme_void()+theme(plot.margin=margin(0,0,0,0,unit="pt"))
  title_grob<-ggplot()+theme_void()+labs(title=title)+
    theme(plot.title=element_text(size=8,hjust=0.5,face="bold"),
          plot.margin=margin(2,2,2,2,unit="pt"))
  narrow_bar<-wrap_plots(list(plot_spacer(),bar,plot_spacer()),nrow=1,
                         widths=c(1,0.6,1))
  wrap_plots(list(title_grob,narrow_bar),ncol=1,heights=c(0.5,1))
}
# legend_cm widened (1.8->2.2) so the two-line gradient title actually fits
# its own column at 8pt (measured ~1.8cm wide, leaving ~0cm clearance before),
# and a third fixed-width spacer column added so any residual sub-mm overflow
# bleeds into blank space instead of clipping against the canvas's true right
# edge — that edge clipping only appeared at large/extreme download canvases,
# never in-app, since the legend column sits flush against the rightmost edge.
gsea_with_legend <- function(main_plot,low_col="#3498db",high_col="#e74c3c",
                              title="-log10(padj)\n(low → high)",legend_cm=2.2,bar_cm=2.5,right_pad_cm=1.2){
  legend_plot<-wrap_plots(list(gsea_legend_bar(low_col,high_col,title),plot_spacer()),ncol=1)+
    plot_layout(heights=grid::unit(c(bar_cm,1),c("cm","null")))
  (main_plot + legend_plot + plot_spacer() +
    plot_layout(widths=grid::unit(c(1,legend_cm,right_pad_cm),c("null","cm","cm"))))
}
# gsea_with_legend()'s legend+pad columns are fixed-cm and get carved out of
# whatever total width is requested, on top of (not part of) barflip_dims'
# estimate for the main plot — so GSEA call sites add this to dd$w to keep
# the main panel from being squeezed further than ORA's (which has no
# external legend column).
gsea_legend_extra_in <- (2.2+1.2)/2.54

# Dot/bar plot of top pathways. Pathway names render horizontally (no rotation
# needed — coord_flip/reorder already puts them on the left, unlike cg_vln/
# cg_dot's rotated bottom labels), so only left-margin + font scaling apply.
# Sizing (barflip_dims) and font floor/ceiling (mkr_font) reuse the same
# helpers as comp_bar/build_mkr_heatmap rather than new constants.
trunc_pw <- function(labels, max_len=100) {
  ifelse(nchar(labels) > max_len, paste0(substr(labels,1,max_len-3),"..."), labels)
}

# Shared msigdbr fetch — called by both gsea_pathways and ora_pathways reactives
# so the if/else collection logic is never duplicated.
fetch_msig <- function(species_str, collection_choice) {
  msig_res <- if(startsWith(collection_choice,"C5_")){
    go_subcat <- paste0("GO:",sub("C5_","",collection_choice))
    msigdbr(species=species_str,collection="C5",subcollection=go_subcat)
  } else if(collection_choice=="H"){
    msigdbr(species=species_str,collection="H")
  } else if(grepl(":",collection_choice)){
    parts<-strsplit(collection_choice,"_")[[1]]; main_coll<-parts[1]
    sub_parts<-strsplit(parts[2],":")[[1]]; sub_coll<-sub_parts[1]
    sub_cat<-if(length(sub_parts)>1) sub_parts[2] else NULL
    if(!is.null(sub_cat)) msigdbr(species=species_str,collection=main_coll,subcollection=paste0(sub_coll,":",sub_cat))
    else msigdbr(species=species_str,collection=main_coll,subcollection=sub_coll)
  } else {
    msigdbr(species=species_str,collection=collection_choice)
  }
  pathways_list <- list(); term2gene <- NULL
  if(!is.null(msig_res) && nrow(msig_res)>0){
    pathways_list <- split(msig_res$gene_symbol, msig_res$gs_name)
    pathways_list <- lapply(pathways_list, function(v) unique(na.omit(v)))
    pathways_list <- pathways_list[vapply(pathways_list,length,1L)>0]
    term2gene <- msig_res %>% dplyr::select(gs_name,gene_symbol) %>%
      dplyr::rename(term=gs_name, gene=gene_symbol)
  }
  list(pathways_list=pathways_list, term2gene=term2gene)
}

# ORA bar/dot plot — single combined view (both directions on one canvas).
# Diverging layout mirrors pe_build_bar's NES>0 fill convention: terms are
# given a signed -log10(p.adjust) (positive = higher in s1, negative = higher
# in s2) so reorder() naturally groups+orders by direction without needing
# two panels. Mirrors gsea_make_plot conventions otherwise (trunc_pw, mkr_font,
# barflip_dims sizing, axis.text.y bold floor 9, min_text_theme, pe color codes).
ora_make_plot <- function(combined_df, comp_label, s1, s2, plot_type="bar") {
  n    <- nrow(combined_df)
  disp <- trunc_pw(combined_df$Description)
  combined_df$term_disp <- disp
  fsz  <- mkr_font(n, 9, 11, 200)
  lab_up <- paste0("Higher in ", s1); lab_dn <- paste0("Higher in ", s2)
  combined_df$signed <- ifelse(combined_df$direction==lab_up,
                                -log10(combined_df$p.adjust),
                                 log10(combined_df$p.adjust))
  cols <- setNames(c("#e74c3c","#3498db"), c(lab_up,lab_dn))
  title_str <- paste0("ORA — ", comp_label)
  base <- if(identical(plot_type,"dot")){
    ggplot(combined_df, aes(x=signed, y=reorder(term_disp,signed),
                            size=Count, color=direction)) +
      geom_point(alpha=0.85) +
      scale_color_manual(values=cols, name=NULL) +
      labs(title=title_str, x="-log10(adj p-value)  [sign = direction]", y=NULL, size="Gene count")
  } else {
    ggplot(combined_df, aes(x=reorder(term_disp,signed), y=signed, fill=direction)) +
      geom_hline(yintercept=0, color="black", linewidth=0.8) +
      geom_col(alpha=0.85) +
      coord_flip() +
      scale_y_continuous(expand=expansion(mult=c(0.02,0.08))) +
      scale_fill_manual(values=cols, name=NULL) +
      labs(title=title_str, x=NULL, y="-log10(adj p-value)  [sign = direction]")
  }
  base + theme_minimal(base_size=12) +
    theme(axis.text.y = element_text(size = fsz, face = "bold"),
          plot.margin=margin(t=10,r=10,b=10,l=max(80,max_label_len(disp)*5),unit="pt"),
          # "panel" (ggplot2 default) anchors the title to the panel's own
          # left edge, which can sit most of the way across the canvas once
          # long category labels push the panel far to the right — at that
          # point a normal-length title has nowhere to go but past the
          # canvas edge. "plot" anchors it to the full figure width instead.
          plot.title.position="plot") +
    min_text_theme
}

gsea_make_plot <- function(top_df,plot_type,s1,s2,comp_label){
  subt<-paste(strwrap(paste0("Positive NES = higher in ",s1,"  |  Negative NES = higher in ",s2),width=70),collapse="\n")
  labels<-top_df$pathway; n<-nrow(top_df)
  disp<-trunc_pw(labels)
  top_df$pathway_disp<-disp
  fsz<-mkr_font(n,9,11,200)
  base<-if(identical(plot_type,"bar")){
    ggplot(top_df,aes(x=reorder(pathway_disp,NES),y=NES,fill=-log10(padj)))+
      geom_hline(yintercept=0,color="black",linewidth=0.8)+
      geom_col(alpha=0.85)+coord_flip()+
      scale_y_continuous(expand=expansion(mult=c(0.02,0.08)))+
      scale_fill_gradient(low="#3498db",high="#e74c3c",guide="none")+
      labs(title=comp_label,subtitle=subt,x=NULL,y="NES")
  } else {
    ggplot(top_df,aes(x=NES,y=reorder(pathway_disp,NES),size=size,color=-log10(padj)))+
      geom_point(alpha=0.85)+
      scale_color_gradient(low="#3498db",high="#e74c3c",guide="none")+
      labs(title=comp_label,subtitle=subt,x="NES",y=NULL,size="Gene set size")
  }
  base+theme_minimal(base_size=12)+
    theme(axis.text.y = element_text(size = fsz, face = "bold"),
          plot.margin=margin(t=10,r=10,b=10,l=max(80,max_label_len(disp)*5),unit="pt"),
          plot.subtitle=element_text(size=10,color="#555"),
          # see ora_make_plot's comment — same panel-vs-plot title anchor fix
          plot.title.position="plot")+
    min_text_theme
}

theme_pub_volcano <- function(base_size=13) {
  theme_minimal(base_size=base_size)+theme(
    axis.text=element_text(size=12,color="black"),
    axis.title=element_text(size=13,face="bold",color="black"),
    axis.line=element_line(color="black",linewidth=0.6),
    panel.grid.minor=element_blank(),legend.position="top",
    plot.title=element_text(size=14,face="bold",hjust=0.5))
}

make_volcano <- function(de_df,lfc_thresh=1,pval_thresh=0.05,top_n=10,title="DE Volcano",s1="A",s2="B") {
  df<-as.data.frame(de_df)
  if(!"gene"%in%colnames(df)) df$gene<-rownames(df)
  pval_col<-if("p_val_adj"%in%colnames(df)) "p_val_adj" else "p_val"
  lfc_col<-"avg_log2FC"
  df<-df[is.finite(df[[lfc_col]])&is.finite(df[[pval_col]])&df[[pval_col]]>0,,drop=FALSE]
  if(nrow(df)==0){plot.new();text(0.5,0.5,"No finite values to plot");return()}
  df$yval<- -log10(df[[pval_col]])
  df$sig <- df[[pval_col]]<=pval_thresh & abs(df[[lfc_col]])>=lfc_thresh
  df$dir <- ifelse(!df$sig,"NS",ifelse(df[[lfc_col]]>0,"Up","Down"))
  up_df<-df[df$dir=="Up",,drop=FALSE]; up_df<-up_df[order(-up_df[[lfc_col]]),,drop=FALSE]
  dn_df<-df[df$dir=="Down",,drop=FALSE]; dn_df<-dn_df[order(dn_df[[lfc_col]]),,drop=FALSE]
  lab_df<-rbind(head(up_df,top_n),head(dn_df,top_n))
  ggplot(df,aes(x=.data[[lfc_col]],y=yval,color=dir))+
    geom_point(alpha=0.55,size=1.8)+
    scale_color_manual(values=c("Up"="#e74c3c","Down"="#3498db","NS"="#aaaaaa"),
                       labels=c("Up"=paste0("Higher in ",s1),"Down"=paste0("Higher in ",s2),"NS"="NS"))+
    {if(nrow(lab_df)>0) geom_text_repel(data=lab_df,aes(label=gene),size=3,color="black",max.overlaps=40,box.padding=0.4)}+
    geom_hline(yintercept=-log10(pval_thresh),linetype="dashed",color="black",linewidth=0.6)+
    geom_vline(xintercept=c(-lfc_thresh,lfc_thresh),linetype="dashed",color="black",linewidth=0.6)+
    labs(title=title,x=paste0("avg log\u2082FC  [log\u2082(",s1," / ",s2,")]"),
         y="-log\u2081\u2080(adj p-value)",color=NULL)+
    theme_pub_volcano()
}

# =============================================================================
# UI
# =============================================================================
ui <- fluidPage(
  titlePanel("scRNA-seq Analysis — SCGUI"),
  tags$head(tags$style(HTML("
    .sh{color:#2c3e50;border-bottom:2px solid #3498db;padding-bottom:4px;
        margin-top:14px;margin-bottom:10px;font-size:16px;font-weight:600}
    .ib{background:#eaf4fb;border-left:4px solid #3498db;padding:7px 11px;
        border-radius:4px;margin-bottom:8px;font-size:12px}
    .wb{background:#fef9e7;border-left:4px solid #f39c12;padding:7px 11px;
        border-radius:4px;margin-bottom:8px;font-size:12px}
    .ob{background:#eafaf1;border-left:4px solid #27ae60;padding:7px 11px;
        border-radius:4px;margin-bottom:8px;font-size:12px}
    .chip{display:inline-block;background:#d6eaf8;border:1px solid #3498db;
          border-radius:10px;padding:1px 9px;margin:2px;font-size:11px}
    hr{border-top:1px solid #ddd}
    .well{padding:10px;margin-bottom:8px}
  "))),

  sidebarLayout(
    # =========================================================================
    # SIDEBAR — Load + Export only
    # =========================================================================
    sidebarPanel(width=3,
      h4("1. Load Data",class="sh"),
      selectInput("species","Species:",c("Human"="human","Mouse"="mouse","Other"="other")),
      radioButtons("data_format","Format:",
                   c("10X MTX"="mtx","10X H5"="h5","Directory"="directory")),
      checkboxInput("multi_dataset","Multiple Datasets",FALSE),
      # single
      conditionalPanel("!input.multi_dataset&&input.data_format=='mtx'",
        fileInput("bc1","Barcodes:",accept=c(".tsv",".txt",".gz")),
        fileInput("gn1","Genes:",   accept=c(".tsv",".txt",".gz")),
        fileInput("mx1","Matrix:",  accept=c(".mtx",".gz"))
      ),
      conditionalPanel("!input.multi_dataset&&input.data_format=='h5'",
        fileInput("h5_single","H5 file:",accept=c(".h5",".h5ad"))
      ),
      conditionalPanel("!input.multi_dataset&&input.data_format=='directory'",
        textInput("data_dir","Dir path:",placeholder="/path/to/filtered_feature_bc_matrix/")
      ),
      # multi sequential queue
      conditionalPanel("input.multi_dataset",
        div(class="ib","Add samples one at a time, then Load."),
        uiOutput("seq_ui"), uiOutput("chips_ui"), hr()
      ),
      conditionalPanel("output.data_loaded",
        h4("Active Samples",class="sh"),
        uiOutput("active_samples_ui"),
        div(class="ib","Controls which samples are shown in UMAP, violin, FeaturePlot, and Comparative DE DotPlot views below. Does not affect exports, statistics, or composition."),
        hr()
      ),
      actionButton("load_data","Load Data",class="btn-primary btn-block"),
      hr(),
      fileInput("rds_file","Load saved .rds object",accept=".rds"),
      actionButton("load_rds","Load RDS",class="btn-primary"),
      textOutput("rds_status"),
      hr(),
      conditionalPanel("output.data_loaded",
        h4("Export",class="sh"),
        downloadButton("dl_seurat","Seurat Object (.rds)"),br(),br(),
        downloadButton("dl_meta_csv","Full Metadata (.csv)")
      )
    ),

    # =========================================================================
    # MAIN PANEL
    # =========================================================================
    mainPanel(width=9, tabsetPanel(id="tabs",

      # -----------------------------------------------------------------------
      # 1 DATA SUMMARY
      # -----------------------------------------------------------------------
      tabPanel("Data Summary",
        h3("Dataset Overview",class="sh"),
        uiOutput("cond_ui"),
        verbatimTextOutput("data_summary"),
        conditionalPanel("output.is_multi",
          h4("Cells per Sample"),
          fluidRow(column(9,plotOutput("samp_dist",height="270px")),
                   column(3,br(),downloadButton("dl_samp_dist","Download")))
        ),
        h4("Cell Metadata — all cells (paginated)"),
        DTOutput("meta_dt"),br(),
        downloadButton("dl_meta_csv2","Download Metadata (.csv)")
      ),

      # -----------------------------------------------------------------------
      # 2 COUNTS MATRIX
      # -----------------------------------------------------------------------
      tabPanel("Counts Matrix",
        conditionalPanel("output.data_loaded",
          h3("Expression Matrix",class="sh"),
          div(class="ib",
            "Showing the joined matrix across all samples. ",
            "Normalized = log-normalized counts [log(counts/total × 10,000 + 1)]. ",
            "Download exports as TSV to handle large cell counts."
          ),
          fluidRow(
            column(3,selectInput("mat_type","Type:",
                     c("Raw Counts"="counts","Normalized"="data","Scaled"="scale.data"))),
            column(3,numericInput("mat_genes","Genes to show:",100,min=10,max=5000)),
            column(3,numericInput("mat_cells","Cells to show:", 100,min=10,max=5000)),
            column(3,uiOutput("mat_samp_ui"))
          ),
          actionButton("update_mat","Update View",class="btn-primary"),br(),br(),
          DTOutput("mat_dt"),br(),
          downloadButton("dl_mat","Download Full Matrix (.tsv)")
        )
      ),

      # -----------------------------------------------------------------------
      # 3 QC METRICS — controls live here
      # -----------------------------------------------------------------------
      tabPanel("QC Metrics",
        h3("Quality Control",class="sh"),
        conditionalPanel("output.obj_raw_avail",
          wellPanel(
            h5("Filter Parameters"),
            numericInput("min_feat","Min Features:",value=200,min=0,max=10000),
            numericInput("max_feat","Max Features:",value=6000,min=0,max=15000),
            numericInput("max_mt","Max % MT:",value=5,min=0,max=100),
            actionButton("apply_qc","Apply QC Filters",class="btn-success btn-block")
          )
        ),
        conditionalPanel("output.qc_applied",
          radioButtons("qc_view","Show:",c("Before"="before","After"="after"),
                       selected="before",inline=TRUE)
        ),
        fluidRow(column(9,plotOutput("qc_vln",height="420px")),
                 column(3,br(),br(),downloadButton("dl_qc_vln","Download"))),
        fluidRow(column(9,plotOutput("qc_scat",height="420px")),
                 column(3,br(),br(),downloadButton("dl_qc_scat","Download"))),
        conditionalPanel("output.is_multi",
          h4("QC by Sample"),
          fluidRow(column(9,plotOutput("qc_samp",height="420px")),
                   column(3,br(),br(),downloadButton("dl_qc_samp","Download")))
        ),
        conditionalPanel("output.qc_applied",
          h4("Filter Summary"),verbatimTextOutput("filt_summ"))
      ),

      # -----------------------------------------------------------------------
      # 4 PCA
      # -----------------------------------------------------------------------
      tabPanel("PCA Analysis",
        conditionalPanel("output.qc_applied",
          h3("PCA",class="sh"),
          div(class="ib",
            "PCA computed on merged data before integration. ",
            "Sample separation here is expected — use the Elbow Plot to choose Dims. ",
            "Heatmaps show the top genes driving each PC — verify PCs capture real biology before integration."
          ),
          wellPanel(fluidRow(
            column(4,numericInput("n_load","PCs for loadings:",4,min=1,max=50)),
            column(4,numericInput("n_heat","PCs for heatmaps:",12,min=1,max=50))
          )),
          h4("PCA plot"),
          fluidRow(column(9,plotOutput("pca_p",height="460px")),column(3,br(),downloadButton("dl_pca_p","Download"))),
          conditionalPanel("output.is_multi",
            h4("PCA by sample"),
            fluidRow(column(9,plotOutput("pca_samp",height="460px")),column(3,br(),downloadButton("dl_pca_samp","Download")))
          ),
          h4("Elbow Plot — choose Dims here"),
          fluidRow(column(9,plotOutput("elbow",height="400px")),column(3,br(),downloadButton("dl_elbow","Download"))),
          h4("Feature Loadings"),
          fluidRow(column(9,plotOutput("pca_load",height="520px")),column(3,br(),downloadButton("dl_pca_load","Download"))),
          h4("PC Heatmaps"),
          fluidRow(column(9,plotOutput("pca_heat",height="1100px")),column(3,br(),downloadButton("dl_pca_heat","Download")))
        )
      ),

      # -----------------------------------------------------------------------
      # 5 INTEGRATION — controls live here
      # -----------------------------------------------------------------------
      tabPanel("Integration",
        h3("Integration",class="sh"),
        conditionalPanel("output.qc_applied",
          wellPanel(
            h5("Integration Method"),
            radioButtons("integ_method","",
                         c("None (PCA)"="none","Harmony"="harmony","CCA"="cca","SCTransform+CCA"="sct_cca"),
                         selected="none",inline=TRUE),
            conditionalPanel("input.integ_method=='harmony'",
              numericInput("harm_theta","Harmony theta:",2,min=0,max=10,step=0.5)
            ),
            conditionalPanel("input.integ_method=='sct_cca'",
              div(class="wb","SCTransform replaces standard normalization. Re-runs from the filtered object.")
            ),
            conditionalPanel("!output.is_multi&&input.integ_method!='none'",
              div(class="wb","Single sample — integration requires 2+ samples. Will use PCA only.")
            ),
            actionButton("run_integ","Run Integration",class="btn-info btn-block")
          )
        ),
        conditionalPanel("!output.qc_applied",
          div(class="wb","Run QC Filtering first (QC Metrics tab).")
        ),
        conditionalPanel("output.integ_done",
          verbatimTextOutput("integ_summ"),
          conditionalPanel("output.is_multi",
            div(class="ib",
              "Both plots below use PCA coordinates — pre-integration. ",
              "Sample separation here is the batch effect being corrected. ",
              "Post-integration UMAP appears in the Clustering tab after you run clustering."
            ),
            conditionalPanel("output.has_integ_preview",
              radioButtons("integ_view","Show:",c("Pre-integration"="pre","Post-integration"="post"),selected="pre",inline=TRUE)
            ),
            fluidRow(
              column(12, h4("Colored by sample"), uiOutput("pre_color_ui"),br(),downloadButton("dl_pre_color","Download"))
            ),
            fluidRow(
              column(12, h4("Split by sample"),   selectInput("pre_split_sample","Sample:",choices=NULL), uiOutput("pre_split_ui"),br(),downloadButton("dl_pre_split","Download"))
            )
          )
        )
      ),

      # -----------------------------------------------------------------------
      # 6 CLUSTERING — controls live here
      # -----------------------------------------------------------------------
      tabPanel("Clustering",
        h3("Clustering",class="sh"),
        conditionalPanel("output.integ_done",
          wellPanel(
            h5("Clustering Parameters"),
            div(class="ib",
              tags$b("Dims")," — PCs fed into graph-building and UMAP. Check the Elbow Plot in the PCA tab.",tags$br(),
              tags$b("Resolution")," — higher = more clusters. Start at 0.5 and adjust based on biology."
            ),
            fluidRow(
              column(4,numericInput("n_dims","Dims:",value=15,min=5,max=50)),
              column(4,numericInput("resolution","Resolution:",value=0.5,min=0.1,max=3,step=0.1))
            ),
            actionButton("run_cluster","Run Clustering",class="btn-info btn-block")
          )
        ),
        conditionalPanel("!output.integ_done",
          div(class="wb","Run Integration first (Integration tab).")
        ),
        conditionalPanel("output.clust_done",
          div(class="ib",uiOutput("clust_reduc_label")),
          fluidRow(
            column(4,uiOutput("clust_ident_ui"),
                   div(style="font-size:11px;color:#555",
                     "After running ScType or applying manual labels in the Annotation tab, this dropdown adds those options.")),
            column(4,radioButtons("umap_col","Color by:",c("Identity"="identity","Sample"="sample"),selected="identity",inline=TRUE))
          ),
          fluidRow(column(9,uiOutput("umap_p_ui")),column(3,br(),downloadButton("dl_umap","Download"))),
          h3("Cluster Statistics"),
          DTOutput("clust_stats_dt"),br(),
          downloadButton("dl_clust_stats","Download (.csv)"),

          hr(),
          h3("Subclustering",class="sh"),
          wellPanel(
            div(class="ib","Subcluster the cluster(s) currently selected from the active identity — works on cluster numbers, ScType labels, manual labels, or a prior subcluster round."),
            uiOutput("sub_target_clusters_ui"),
            actionButton("prep_subcluster","1. Prep Subset & Run PCA",class="btn-warning"),
            conditionalPanel("output.subprep_done",
              hr(),
              radioButtons("sub_integ_method","Subcluster Integration:",
                           c("Match global","None","Harmony","CCA","SCTransform+CCA"),
                           selected="Match global"),
              conditionalPanel("input.sub_integ_method=='Harmony'",
                numericInput("sub_harm_theta","Subcluster Harmony theta:",value=2,min=0,max=10,step=0.5)
              ),
              actionButton("integ_subcluster","2. Integrate",class="btn-info")
            ),
            conditionalPanel("output.sub_integ_done",
              hr(),
              numericInput("sub_dims","Subcluster Dims:",value=10,min=2,max=30),
              numericInput("sub_resolution","Subcluster Resolution:",value=0.5,min=0.1,max=3,step=0.1),
              actionButton("run_subcluster","3. Cluster",class="btn-info")
            )
          ),
          conditionalPanel("output.subprep_done",
            uiOutput("sub_diag_header"),
            fluidRow(
              column(6,h4("Elbow Plot"),uiOutput("sub_elbow_ui"),br(),downloadButton("dl_sub_elbow","Download")),
              column(6,h4("PC Heatmap"),uiOutput("sub_pca_heat_ui"),br(),downloadButton("dl_sub_pca_heat","Download"))
            ),
            conditionalPanel("output.is_multi",
              h4("Unintegrated UMAP by Sample"),
              conditionalPanel("output.sub_integ_done",
                radioButtons("sub_integ_view","Show:",c("Pre-integration"="pre","Post-integration"="post"),selected="pre",inline=TRUE)
              ),
              uiOutput("sub_umap_pre_ui"),br(),downloadButton("dl_sub_umap_pre","Download"),
              div(class="ib","If cells separate cleanly by sample above with no biological interpretation, consider reintegrating below. If they mix reasonably or the separation looks biologically meaningful (e.g. activation state), you likely don't need to.")
            )
          ),
          conditionalPanel("output.subclust_done",
            hr(),
            h4("Subcluster Result"),
            fluidRow(column(9,uiOutput("sub_umap_ui")),column(3,br(),downloadButton("dl_sub_umap","Download"))),
            h5("Subcluster Sizes"),
            DTOutput("sub_clust_dt"),
            downloadButton("dl_sub_clust_dt","Download (.csv)"),
            br(),
            h5("Subcluster Sizes by Sample"),
            DTOutput("sub_clust_sample_dt"),
            downloadButton("dl_sub_clust_sample_dt","Download (.csv)"),
            br(),
            conditionalPanel("output.sub_multi_parent",
              h5("Subcluster Composition by Parent Origin (multi-cluster subclustering only)"),
              DTOutput("sub_clust_origin_dt"),
              downloadButton("dl_sub_clust_origin_dt","Download (.csv)"),
              br()
            ),
            hr(),
            h4("Subcluster Top Markers"),
            div(class="ib","Compares subclusters against EACH OTHER (within this subset), not against the whole dataset — use this to find what distinguishes subclusters that all washed out as 'the same thing' when compared to the rest of the data."),
            wellPanel(
              fluidRow(
                column(3,numericInput("sub_top_n","Top N per cluster:",10,min=1,max=500)),
                column(3,numericInput("sub_min_pct","Min pct expressed:",0.25,min=0.05,max=1,step=0.05)),
                column(3,numericInput("sub_min_lfc","Min log2FC:",0.25,min=0.05,max=5,step=0.05)),
                column(3,br(),actionButton("find_sub_markers","Find Subcluster Markers",class="btn-primary btn-block"))
              )
            ),
            conditionalPanel("output.sub_markers_done",
              tabsetPanel(
                tabPanel("Heatmap",
                  fluidRow(column(9,uiOutput("sub_mkr_heat_ui")),column(3,br(),downloadButton("dl_sub_mkr_heat","Download")))
                ),
                tabPanel("Top Markers Table",
                  DTOutput("sub_mkr_top_dt"),br(),
                  fluidRow(column(4,downloadButton("dl_sub_mkr_top","Top N (.csv)")),
                           column(4,downloadButton("dl_sub_mkr_all","All Markers (.csv)")))
                )
              )
            ),
            actionButton("merge_subcluster","Merge into main object",class="btn-success"),
            actionButton("discard_subcluster","Discard and retry",class="btn-danger")
          )
        )
      ),

      # -----------------------------------------------------------------------
      # 7 CLUSTER QC
      # -----------------------------------------------------------------------
      tabPanel("Cluster QC",
        conditionalPanel("output.clust_done",
          h3("QC per Cluster",class="sh"),
          fluidRow(column(9,plotOutput("cqc_vln",height="460px")),column(3,br(),br(),downloadButton("dl_cqc_vln","Download"))),
          fluidRow(column(9,plotOutput("cqc_scat",height="460px")),column(3,br(),br(),downloadButton("dl_cqc_scat","Download")))
        ),
        conditionalPanel("!output.clust_done", div(class="wb","Run Clustering first."))
      ),

      # -----------------------------------------------------------------------
      # 8 GENE EXPRESSION
      # -----------------------------------------------------------------------
      tabPanel("Gene Expression",
        conditionalPanel("output.clust_done",
          h3("Gene Expression",class="sh"),
          div(class="ib",
            "Exploratory view on the full integrated dataset. ",
            "Y-axis = log-normalized expression: log(counts/total × 10,000 + 1). ",
            "UMAP coordinates reflect your chosen integration method. ",
            "For condition comparisons (ctrl vs lps) → Comparative Analysis tab."
          ),
          wellPanel(fluidRow(
            column(4,textInput("gene_in","Gene(s) — comma-separated:","")),
            column(3,selectInput("gene_type","Plot type:",c("Feature Plot"="feature","Violin Plot"="violin","Dot Plot"="dot"))),
            column(3,uiOutput("gene_ident_ui")),
            column(2,br(),actionButton("plot_gene","Plot",class="btn-warning btn-block"))
          ),
          fluidRow(
            column(6,uiOutput("gene_clust_filter_ui"))
          )),
          fluidRow(column(9,div(style="overflow-x:auto;",uiOutput("gene_p_ui"))),column(3,br(),downloadButton("dl_gene","Download")))
        ),
        conditionalPanel("!output.clust_done", div(class="wb","Run Clustering first."))
      ),

      # -----------------------------------------------------------------------
      # 9 CLUSTER MARKERS
      # -----------------------------------------------------------------------
      tabPanel("Cluster Markers",
        conditionalPanel("output.clust_done",
          h3("Cluster Markers",class="sh"),
          wellPanel(
            fluidRow(
              column(3,numericInput("top_n","Top N per cluster:",10,min=1,max=500)),
              column(3,numericInput("min_pct","Min pct expressed:",0.25,min=0.05,max=1,step=0.05)),
              column(3,numericInput("min_lfc","Min log2FC:",0.25,min=0.05,max=5,step=0.05)),
              column(3,br(),actionButton("find_markers","Find All Markers",class="btn-primary btn-block"))
            ),
            textOutput("mkr_group_by_txt"),
            div(class="ib",
              tags$b("Min pct")," = fraction of cluster cells expressing the gene (0.25 = 25%). ",
              tags$b("Min log2FC")," = fold-change threshold (0.25 ≈ 1.19×). Seurat tutorial defaults."
            )
          ),
          conditionalPanel("output.markers_done",
            tabsetPanel(
              tabPanel("Heatmap",
                fluidRow(column(9,uiOutput("mkr_heat_ui")),column(3,br(),downloadButton("dl_mkr_heat","Download")))
              ),
              tabPanel("Top Markers Table",
                DTOutput("mkr_top_dt"),br(),
                fluidRow(column(4,downloadButton("dl_mkr_top","Top N (.csv)")),
                         column(4,downloadButton("dl_mkr_all","All Markers (.csv)")))
              ),
              tabPanel("All Genes by Cluster",
                div(class="ib","All markers for a selected cluster, sorted by log2FC. Set N=0 for all."),
                fluidRow(column(4,uiOutput("mkr_clust_sel_ui")),column(4,numericInput("mkr_show_n","Show top N (0=all):",0,min=0))),
                DTOutput("mkr_cluster_dt"),br(),
                downloadButton("dl_mkr_cluster","Download this cluster (.csv)")
              )
            )
          )
        ),
        conditionalPanel("!output.clust_done", div(class="wb","Run Clustering first."))
      ),

      # -----------------------------------------------------------------------
      # 10 ANNOTATION
      # -----------------------------------------------------------------------
      tabPanel("Annotation",
        conditionalPanel("output.clust_done",
          h3("Cell Type Annotation",class="sh"),
          tabsetPanel(

            tabPanel("ScType (Automated)",
              h4("ScType — Database-Driven Annotation"),
              div(class="ib",
                "ScType matches your cluster expression profiles against a curated cell type marker database. ",
                "The table shows which genes from the database were expressed in each cluster and drove the annotation. ",
                "Higher score = stronger match. Clusters below 25% of max score are labelled Unknown."
              ),
              wellPanel(
                selectInput("tissue","Tissue / Organ:",
                  c("Immune system","Brain","Lung","Liver","Pancreas","Kidney","Eye","Heart",
                    "Adrenal","Intestine","Muscle","Placenta","Spleen","Stomach","Thymus"),
                  selected="Immune system"),
                actionButton("run_sctype","Run ScType",class="btn-success")
              ),
              conditionalPanel("output.sctype_done",
                fluidRow(column(9,uiOutput("sctype_umap_ui")),column(3,br(),downloadButton("dl_sctype_umap","Download"))),
                h5("Results — Cluster assignments and marker genes used"),
                DTOutput("sctype_tbl"),br(),
                downloadButton("dl_sctype_tbl","Download (.csv)")
              )
            ),

            tabPanel("Signature Scoring",
              div(class="ib",
                "Define custom gene sets from the literature (CellMarker, PanglaoDB, etc.). ",
                "AddModuleScore scores each cluster against each signature. Results guide manual annotation — ",
                "top-scoring signature per cluster appears in the Manual Annotation tab."
              ),
              wellPanel(
                h5("Upload Signatures CSV"),
                div(class="ib", "CSV must contain columns named Signature and Genes. Genes should be comma-separated."),
                fileInput("sig_csv", "Signature CSV:", accept = c(".csv")),
                actionButton("load_sig_csv", "Load Signatures CSV", class = "btn-primary btn-sm"),
                hr(),
                h5("Add a Signature Manually"),
                fluidRow(
                  column(4, textInput("sig_name","Signature name:","")),
                  column(6, textInput("sig_genes","Genes (comma-separated):",""))
                ),
                fluidRow(
                  column(4, actionButton("add_sig","Add Signature",class="btn-info btn-sm")),
                  column(4, actionButton("clear_sigs","Clear All",class="btn-danger btn-sm"))
                )
              ),
              uiOutput("sig_list_ui"),
              br(),
              actionButton("run_sig","Score All Clusters",class="btn-success"),
              textOutput("sig_group_by_txt"),
              br(),br(),
              conditionalPanel("output.sig_done",
                fluidRow(
                  column(9, uiOutput("sig_heat_ui")),
                  column(3, br(), downloadButton("dl_sig_heat","Download PNG"), br(),br(),
                              downloadButton("dl_sig_csv","Download CSV"))
                ),
                h5("Raw Average Scores per Cluster"),
                DTOutput("sig_tbl"),
                hr(),
                h5("Signature Gene Inspector"),
                uiOutput("sig_view_sel_ui"),
                DTOutput("sig_gene_dt"),
                br(),
                downloadButton("dl_sig_gene_csv","Download Gene Contributions — All Signatures (.csv)")
              )
            ),

            tabPanel("Manual Annotation",
              div(class="ib",
                "Type a free-text label for each cluster. ScType suggestion and top 5 marker genes are shown for guidance. ",
                "Marker genes are only available after running 'Find All Markers'. ",
                "After applying, 'Manual Label' becomes available in all identity dropdowns throughout the app."
              ),
              uiOutput("manual_annot_ui"),
              br(),
              fluidRow(
                column(4,actionButton("apply_manual_annot","Apply Labels",class="btn-success")),
                column(4,actionButton("reset_manual_annot","Reset",class="btn-warning")),
                column(4,actionButton("refresh_mannot","Refresh Hints",class="btn-default btn-sm",
                                      title="Reload hints after running ScType, Find Markers, or Signature Scoring"))
              ),
              textOutput("mannot_group_by_txt"),
              verbatimTextOutput("manual_annot_status"),
              conditionalPanel("output.manual_annot_done",
                h4("Manual Annotation UMAP"),
                fluidRow(column(9,uiOutput("manual_annot_umap_ui")),column(3,br(),downloadButton("dl_manual_annot_umap","Download")))
              )
            )
          )
        ),
        conditionalPanel("!output.clust_done", div(class="wb","Run Clustering first."))
      ),

      # -----------------------------------------------------------------------
      # 11 COMPARATIVE ANALYSIS
      # -----------------------------------------------------------------------
      tabPanel("Comparative Analysis",
        conditionalPanel("output.clust_done",
          conditionalPanel("!output.is_multi", div(class="wb","Comparative Analysis requires 2+ samples.")),
          conditionalPanel("output.is_multi",
            h3("Comparative Analysis",class="sh"),
            fluidRow(
              column(3,selectInput("comp_ident","Color / group by:",c("Cluster Number"="cluster"),selected="cluster")),
              column(9,helpText("Identity choice applies to all plots in this tab."))
            ),
            tabsetPanel(

              tabPanel("UMAP",
                h4("UMAP Overview",class="sh"),
                tabsetPanel(
                  tabPanel("Colored by Sample",fluidRow(column(9,uiOutput("cu_samp_ui")),column(3,br(),downloadButton("dl_cu_samp","Download")))),
                  tabPanel("Split by Sample",  fluidRow(column(9,selectInput("cu_split_sample","Sample:",choices=NULL),
                           radioButtons("cu_split_color_by","Color by:",c("Identity"="identity","Sample"="sample"),selected="identity",inline=TRUE),
                           uiOutput("cu_split_ui")),column(3,br(),downloadButton("dl_cu_split","Download")))),
                  tabPanel("Colored by Identity",fluidRow(column(9,uiOutput("cu_ident_ui")),column(3,br(),downloadButton("dl_cu_ident","Download"))))
                )
              ),

              tabPanel("Composition",
                h4("Cluster Composition by Sample",class="sh"),
                div(class="ib","Proportion of cells from each sample per cluster or cell type."),
                fluidRow(column(9,uiOutput("comp_bar_ui")),column(3,br(),br(),downloadButton("dl_comp_bar","Download"))),
                br(),downloadButton("dl_comp_csv","Download (.csv)"),
                hr(),
                h4("Per-Sample Cell Type Proportions"),
                div(class="ib","For each sample: what percentage of that sample's cells belong to each cell type."),
                fluidRow(column(9,div(style="overflow-x:auto;",uiOutput("comp_bar_persample_ui"))),column(3,br(),br(),downloadButton("dl_comp_bar_persample","Download PNG"),br(),br(),downloadButton("dl_comp_persample_csv","Download CSV")))
              ),

              tabPanel("Gene Expression",
                h4("Comparative Gene Expression",class="sh"),
                div(class="ib",uiOutput("comp_reduc_label")),
                div(class="ib","Y-axis = log-normalized expression. Split plots show each sample side by side."),
                wellPanel(fluidRow(
                  column(8,textInput("comp_genes","Gene(s) — comma-separated:","")),
                  column(4,br(),actionButton("plot_comp","Plot",class="btn-warning btn-block"))
                )),
                tabsetPanel(
                  tabPanel("Split Feature Plot",fluidRow(column(9,uiOutput("cg_feat_ui")),column(3,br(),downloadButton("dl_cg_feat","Download")))),
                  tabPanel("Split Violin Plot", fluidRow(column(9,uiOutput("cg_vln_ui")),
                           column(3,br(),downloadButton("dl_cg_vln","Download"),br(),br(),
                                  radioButtons("cg_vln_scope","Scope:",c("Global"="global","Per-cluster"="percluster","Selected"="selected"),selected="global",inline=TRUE),
                                  conditionalPanel("input.cg_vln_scope=='percluster'",uiOutput("cg_vln_clust_ui")),
                                  conditionalPanel("input.cg_vln_scope=='selected'",uiOutput("cg_vln_clusters_sel_ui"))))),
                  tabPanel("Split Dot Plot",    fluidRow(column(9,div(style="overflow-x:auto;",uiOutput("cg_dot_ui"))),
                           column(3,br(),downloadButton("dl_cg_dot","Download"),br(),br(),
                                  radioButtons("cg_dot_scope","Scope:",c("Global"="global","Per-cluster"="percluster","Selected"="selected"),selected="global",inline=TRUE),
                                  conditionalPanel("input.cg_dot_scope=='percluster'",uiOutput("cg_dot_clust_ui")),
                                  conditionalPanel("input.cg_dot_scope=='selected'",uiOutput("cg_dot_clusters_sel_ui")))))
                )
              ),

              tabPanel("Conserved Markers",
                h4("Conserved Markers",class="sh"),
                div(class="ib",tags$b("FindConservedMarkers")," finds genes that mark a cluster robustly regardless of condition — uses all samples simultaneously."),
                fluidRow(column(4,uiOutput("cons_clust_ui")),column(4,br(),actionButton("run_cons","Run",class="btn-primary"))),br(),
                conditionalPanel("output.cons_done",DTOutput("cons_dt"),br(),downloadButton("dl_cons","Download (.csv)"))
              ),

              tabPanel("Pseudobulk Scatter",
                h4("Pseudobulk Scatter",class="sh"),
                div(class="ib",
                  "AggregateExpression creates pseudobulk profiles per cell type per sample. ",
                  "Points above the diagonal = higher in Sample A; below = higher in Sample B."
                ),
                fluidRow(column(4,uiOutput("pb_ct_ui")),column(4,textInput("pb_hi","Highlight genes (comma-sep):","")),
                         column(4,br(),actionButton("run_pb","Run Pseudobulk",class="btn-primary"))),br(),
                conditionalPanel("output.pb_done",
                  uiOutput("pb_warn_ui"),
                  fluidRow(column(9,plotOutput("pb_p",height="620px")),column(3,br(),downloadButton("dl_pb","Download")))
                )
              ),

              tabPanel("Condition DE",
                uiOutput("de_samp_ui"),
                h4("Condition DE",class="sh"),
                div(class="ib",uiOutput("de_formula_ui")),
                fluidRow(column(4,uiOutput("de_clust_ui")),column(4,br(),actionButton("run_de","Run DE",class="btn-primary"))),br(),
                div(style="font-size:11px;color:#999","Dev/diagnostic — cached DE comparisons (temporary, remove later):"),
                verbatimTextOutput("de_cache_debug"),
                conditionalPanel("output.de_done",
                  radioButtons("de_view", NULL,
                               choices=c("Results Table"="table", "Volcano Plot"="volcano"),
                               selected="table", inline=TRUE)
                ),
                conditionalPanel("output.de_done && input.de_view == 'table'",
                  DTOutput("de_dt"),br(),downloadButton("dl_de","Download (.csv)")
                ),
                conditionalPanel("output.de_done && input.de_view == 'volcano'",
                  hr(),
                  fluidRow(
                    column(4,numericInput("volc_lfc","Min |log2FC|:",1,min=0,max=10,step=0.1)),
                    column(4,numericInput("volc_pval","Adj. p-value threshold (FDR):",0.05,min=0.001,max=1,step=0.01)),
                    column(4,numericInput("volc_n","Top N labels per side:",10,min=1,max=50))
                  ),
                  fluidRow(column(9,plotOutput("volc_p",height="750px")),column(3,br(),downloadButton("dl_volc","Download")))
                )
              ),

              tabPanel("GSEA",
                h4("Gene Set Enrichment Analysis (GSEA)",class="sh"),
                div(class="ib","Runs fgsea on a cached DE comparison from the Condition DE tab. Species/collection here are independent of the species setting in the sidebar."),
                conditionalPanel("!output.has_de_cache",
                  div(class="wb","No cached DE comparisons yet — run at least one comparison in the Condition DE tab first.")
                ),
                conditionalPanel("output.has_de_cache",
                  wellPanel(
                    fluidRow(
                      column(4,uiOutput("gsea_comparison_ui")),
                      column(4,selectInput("gsea_species","Species:",
                               c("Human"="Homo sapiens","Mouse"="Mus musculus","Rat"="Rattus norvegicus",
                                 "Zebrafish"="Danio rerio","Fly"="Drosophila melanogaster","Worm"="Caenorhabditis elegans"))),
                      column(4,selectInput("gsea_collection","Gene set collection:",
                               c("Hallmark"="H","GO Biological Process"="C5_BP","GO Cellular Component"="C5_CC",
                                 "GO Molecular Function"="C5_MF","KEGG"="C2_CP:KEGG","Reactome"="C2_CP:REACTOME")))
                    ),
                    fluidRow(
                      column(4,radioButtons("gsea_plot_type","Plot type:",c("Dot"="dot","Bar"="bar"),selected="dot",inline=TRUE)),
                      column(4,numericInput("gsea_topn","Top N pathways per direction:",15,min=1,max=50)),
                      column(4,br(),actionButton("run_gsea","Run GSEA",class="btn-primary btn-block"))
                    )
                  ),
                    tabsetPanel(
                    tabPanel("GSEA Results",
                      conditionalPanel("output.gsea_done",
                        fluidRow(
                          column(9,div(style="overflow-x:auto;",uiOutput("gsea_plot_ui"))),
                          column(3,br(),downloadButton("dl_gsea_plot","Download Plot"),br(),br(),
                                    downloadButton("dl_gsea_csv","Download Results (.csv)"))
                        )
                      )
                    ),
                    tabPanel("Pathway Explorer",
                      conditionalPanel("!output.gsea_done",
                        div(class="wb","Run GSEA first (GSEA Results sub-tab) to populate pathway choices.")
                      ),
                      conditionalPanel("output.gsea_done",
                        wellPanel(
                          fluidRow(
                            column(6,selectizeInput("pe_pathways","Select Pathways:",choices=NULL,multiple=TRUE,
                                     options=list(placeholder="Choose pathways...",maxItems=20))),
                            column(6,radioButtons("pe_view","View:",
                                     c("NES Barplot"="barplot","Enrichment Trace"="trace"),
                                     selected="barplot",inline=TRUE))
                          )
                        ),
                        fluidRow(
                          column(9,div(style="overflow-x:auto;",uiOutput("pe_plot_ui"))),
                          column(3,br(),downloadButton("dl_pe_plot","Download"))
                        )
                      )
                    )
                  )
                )
              ),

              tabPanel("ORA",
                h4("Over-Representation Analysis (ORA)",class="sh"),
                div(class="ib","Runs enricher() on significant DE genes from a cached comparison. Species/collection here are independent of all other tabs."),
                conditionalPanel("!output.has_de_cache",
                  div(class="wb","No cached DE comparisons yet — run at least one comparison in the Condition DE tab first.")
                ),
                conditionalPanel("output.has_de_cache",
                  wellPanel(
                    fluidRow(
                      column(4,uiOutput("ora_comparison_ui")),
                      column(4,selectInput("ora_species","Species:",
                               c("Human"="Homo sapiens","Mouse"="Mus musculus","Rat"="Rattus norvegicus",
                                 "Zebrafish"="Danio rerio","Fly"="Drosophila melanogaster","Worm"="Caenorhabditis elegans"))),
                      column(4,selectInput("ora_collection","Gene set collection:",
                               c("Hallmark"="H","GO Biological Process"="C5_BP","GO Cellular Component"="C5_CC",
                                 "GO Molecular Function"="C5_MF","KEGG"="C2_CP:KEGG","Reactome"="C2_CP:REACTOME")))
                    ),
                    fluidRow(
                      column(3,numericInput("ora_pcut","Adj p-value cutoff:",0.05,min=0.001,max=1,step=0.01)),
                      column(3,numericInput("ora_lfc","Min |log2FC|:",0.5,min=0,max=10,step=0.1)),
                      column(3,numericInput("ora_topn","Top N terms:",20,min=1,max=100)),
                      column(3,radioButtons("ora_plot_type","Plot type:",c("Bar"="bar","Dot"="dot"),selected="bar",inline=TRUE))
                    ),
                    fluidRow(
                      column(4,actionButton("run_ora","Run ORA",class="btn-primary btn-block"))
                    ),
                    conditionalPanel("output.ora_has_results",
                      fluidRow(
                        column(4,downloadButton("dl_ora_plot","Download Plot (600 DPI)")),
                        column(4,downloadButton("dl_ora_csv","Download Results (.csv)"))
                      )
                    )
                  ),
                  conditionalPanel("output.ora_done",
                    div(style="overflow-x:auto;",uiOutput("ora_plot_ui"))
                  )
                )
              )
            )
          )
        )
      )
    ))
  )
)

# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  rv <- reactiveValues(
    obj=NULL, obj_raw=NULL, obj_filt=NULL, obj_joined=NULL,
    data_loaded=FALSE, is_multi=FALSE, samp_names=character(0), cond_map=NULL,
    qc_applied=FALSE, n_cells=0, rds_status="",
    integ_done=FALSE, integ_method="none", reduc_use="pca", harm_theta_used=NA,
    clust_done=FALSE,
    markers=NULL, markers_done=FALSE,
    cur_gene_plot=NULL, cur_genes=NULL, cg_state=NULL,
    sctype_done=FALSE, sctype_res=NULL,
    manual_annot_done=FALSE, manual_annot_labels=list(),
    sig_sets=list(), sig_scores=NULL, sig_done=FALSE, sig_gene_table=NULL,
    cons_markers=NULL, cons_done=FALSE,
    pb_obj=NULL, pb_done=FALSE, pb_id_col=NULL,
    de_res=NULL, de_done=FALSE, de_s1="", de_s2="", de_cache=list(),
    gsea_res=NULL, gsea_done=FALSE, gsea_s1="", gsea_s2="", gsea_label="",
    gsea_pl=NULL, gsea_stats=NULL,
    ora_cache=list(), ora_res_up=NULL, ora_res_dn=NULL, ora_done=FALSE,
    ora_s1="", ora_s2="", ora_label="",
    pending=list(),
    sub_obj=NULL, subprep_done=FALSE, sub_integ_done=FALSE, subclust_done=FALSE,
    sub_ru=NULL, sub_was_split=FALSE,
    sub_parent_clusters=NULL, sub_parent_col=NULL,
    sub_markers=NULL, sub_markers_done=FALSE
  )

  # ---- slider sync ----

  # ---- sequential sample queue ----
  output$seq_ui <- renderUI({
    req(input$multi_dataset)
    n<-length(rv$pending)+1; fmt<-input$data_format
    wellPanel(
      h5(paste("Add Sample",n)),
      textInput("new_nm","Sample name:",paste0("Sample_",n)),
      if(fmt=="mtx") tagList(
        fileInput("new_bc","Barcodes:",accept=c(".tsv",".txt",".gz")),
        fileInput("new_gn","Genes:",   accept=c(".tsv",".txt",".gz")),
        fileInput("new_mx","Matrix:",  accept=c(".mtx",".gz"))
      ) else if(fmt=="h5") fileInput("new_h5","H5 file:",accept=c(".h5",".h5ad")),
      fluidRow(
        column(6,actionButton("add_samp","Add Sample",class="btn-success btn-sm")),
        column(6,actionButton("clear_samps","Clear All",class="btn-danger btn-sm"))
      )
    )
  })
  output$chips_ui <- renderUI({
    if(!length(rv$pending)) return(NULL)
    tagList(h6(paste(length(rv$pending),"sample(s) queued:")),
            lapply(rv$pending,function(s) span(class="chip",paste0("✓ ",s$name))))
  })
  output$active_samples_ui <- renderUI({
    req(rv$samp_names)
    checkboxGroupInput("active_samples","Active samples:",choices=rv$samp_names,selected=rv$samp_names)
  })
  observeEvent(input$add_samp,{
    req(input$new_nm); nm<-trimws(input$new_nm)
    if(!nchar(nm)){showNotification("Enter a sample name.",type="warning");return()}
    e<-list(name=nm,fmt=input$data_format)
    if(input$data_format=="mtx"){req(input$new_bc,input$new_gn,input$new_mx);e$bc<-input$new_bc;e$gn<-input$new_gn;e$mx<-input$new_mx}
    else if(input$data_format=="h5"){req(input$new_h5);e$h5<-input$new_h5}
    rv$pending<-c(rv$pending,list(e)); showNotification(paste("Queued:",nm),type="message")
  })
  observeEvent(input$clear_samps,{rv$pending<-list()})

  # ---- LOAD DATA ----
  observeEvent(input$load_data,{
    withProgress(message="Loading data...",value=0,{
      tryCatch({
        if(!input$multi_dataset){
          if(input$data_format=="mtx"){
            req(input$bc1,input$gn1,input$mx1); td<-tempdir()
            file.copy(input$bc1$datapath,file.path(td,"barcodes.tsv"),overwrite=TRUE)
            file.copy(input$gn1$datapath,file.path(td,"genes.tsv"),overwrite=TRUE)
            file.copy(input$mx1$datapath,file.path(td,"matrix.mtx"),overwrite=TRUE)
            dat<-Read10X(data.dir=td,gene.column=2)
          } else if(input$data_format=="h5"){req(input$h5_single);dat<-Read10X_h5(input$h5_single$datapath)
          } else {req(input$data_dir);dat<-Read10X(data.dir=input$data_dir,gene.column=2)}
          obj<-CreateSeuratObject(counts=dat,project="scRNA",min.cells=3,min.features=200)
          obj$sample<-"Sample1"; rv$is_multi<-FALSE
        } else {
          samps<-rv$pending
          if(length(samps)<2){showNotification("Add at least 2 samples.",type="warning");return()}
          sl<-list()
          for(i in seq_along(samps)){
            incProgress(i/length(samps)*0.6,detail=paste("Reading",samps[[i]]$name))
            s<-samps[[i]]
            if(s$fmt=="mtx"){
              td<-file.path(tempdir(),paste0("ds",i));dir.create(td,showWarnings=FALSE)
              file.copy(s$bc$datapath,file.path(td,"barcodes.tsv"),overwrite=TRUE)
              file.copy(s$gn$datapath,file.path(td,"genes.tsv"),overwrite=TRUE)
              file.copy(s$mx$datapath,file.path(td,"matrix.mtx"),overwrite=TRUE)
              d<-Read10X(data.dir=td,gene.column=2)
            } else if(s$fmt=="h5") d<-Read10X_h5(s$h5$datapath)
              else next
            o<-CreateSeuratObject(counts=d,project=s$name,min.cells=3,min.features=200); o$sample<-s$name; sl[[i]]<-o
          }
          incProgress(0.7,detail="Merging...")
          ids<-sapply(sl,function(x) unique(x$sample)[1])
          obj<-merge(sl[[1]],y=sl[-1],add.cell.ids=ids,project="Merged")
          lnames<-names(obj[["RNA"]]@layers)
          if(!any(grepl("\\.",lnames))) obj[["RNA"]]<-split(obj[["RNA"]],f=obj$sample)
          rv$is_multi<-TRUE
        }
        incProgress(0.85,detail="QC metrics...")
        obj[["percent.mt"]]<-PercentageFeatureSet(obj,pattern=get_mt_pat(input$species))
        if(max(obj[["percent.mt"]])==0){
          mg<-grep("^MT-|^mt-|^Mt-",rownames(obj),value=TRUE,ignore.case=TRUE)
          if(length(mg)) obj[["percent.mt"]]<-PercentageFeatureSet(obj,features=mg)
        }
        obj$condition<-obj$sample
        rv$obj<-obj; rv$obj_raw<-obj; rv$obj_joined<-JoinLayers(obj)
        rv$samp_names<-sort(unique(obj$sample)); rv$data_loaded<-TRUE
        rv$qc_applied<-FALSE; rv$integ_done<-FALSE; rv$clust_done<-FALSE
        rv$markers_done<-FALSE; rv$sctype_done<-FALSE; rv$manual_annot_done<-FALSE
        rv$sig_done<-FALSE; rv$sig_sets<-list(); rv$sig_scores<-NULL
        rv$cons_done<-FALSE; rv$pb_done<-FALSE; rv$de_done<-FALSE
        rv$sig_gene_table<-NULL; rv$manual_annot_labels<-list()
        incProgress(1)
        showNotification(paste("Loaded",ncol(obj),"cells from",length(rv$samp_names),"sample(s)"),type="message")
      },error=function(e){showNotification(paste("Load error:",e$message),type="error",duration=15);message(e$message)})
    })
  })

  # ---- LOAD RDS ----
  observeEvent(input$load_rds,{
    req(input$rds_file)
    tryCatch({
      obj <- readRDS(input$rds_file$datapath)
      rv$obj <- obj
      rv$obj_raw <- obj
      rv$obj_filt <- obj
      rv$data_loaded <- TRUE
      rv$qc_applied <- "pca"%in%names(obj@reductions) && "scale.data"%in%SeuratObject::Layers(obj[["RNA"]])
      rv$clust_done <- "umap"%in%names(obj@reductions) && "seurat_clusters"%in%colnames(obj@meta.data)
      rv$sctype_done <- "sctype_classification"%in%colnames(obj@meta.data)
      rv$manual_annot_done <- "manual_annotation"%in%colnames(obj@meta.data)
      rv$manual_annot_labels <- list()
      if(rv$manual_annot_done){
        build_manual_labels <- function(grp_col){
          if(!grp_col %in% colnames(obj@meta.data)) return(NULL)
          grp <- as.character(obj@meta.data[[grp_col]])
          lab <- as.character(obj@meta.data$manual_annotation)
          pairs <- unique(data.frame(grp=grp, lab=lab, stringsAsFactors=FALSE))
          if(anyDuplicated(pairs$grp)) return(NULL)
          setNames(as.list(pairs$lab), pairs$grp)
        }
        labs <- build_manual_labels("seurat_clusters")
        if(is.null(labs)) labs <- build_manual_labels("subcluster_labels")
        if(!is.null(labs)) rv$manual_annot_labels <- labs
      }
      rv$integ_done <- rv$clust_done || any(c("harmony","integrated.cca","integrated.dr") %in% names(obj@reductions))
      rv$reduc_use <- if("harmony"%in%names(obj@reductions)) "harmony"
                       else if("integrated.cca"%in%names(obj@reductions)) "integrated.cca"
                       else if("integrated.dr"%in%names(obj@reductions)) "integrated.dr"
                       else "pca"
      rv$integ_method <- switch(rv$reduc_use, harmony="harmony", integrated.cca="cca",
                                 integrated.dr="sct_cca", "none")
      updateRadioButtons(session, "integ_method", selected=rv$integ_method)
      rv$harm_theta_used <- NA
      rv$samp_names <- nat_sort(unique(as.character(obj@meta.data[["sample"]]%||%obj@meta.data[["orig.ident"]])))
      rv$is_multi <- length(rv$samp_names)>1
      rv$n_cells <- ncol(obj)
      n_clust <- if(rv$clust_done) length(unique(obj$seurat_clusters)) else NA
      rv$rds_status <- paste0(
        "Loaded ",format(rv$n_cells,big.mark=",")," cells | ",
        length(rv$samp_names)," sample(s)",
        if(rv$clust_done) paste0(" | ",n_clust," clusters") else "",
        if(rv$sctype_done) " | ScType ✓" else "",
        if(rv$manual_annot_done) " | Manual annotation ✓" else "")
      showNotification("RDS loaded successfully.",type="message")
    },error=function(e){showNotification(paste("RDS load error:",e$message),type="error",duration=15);message(e$message)})
  })
  output$rds_status <- renderText(rv$rds_status%||%"")

  # ---- CONDITION UI ----
  output$cond_ui <- renderUI({
    req(rv$data_loaded,rv$is_multi)
    tagList(
      h4("Optional: Assign Condition Groups"),
      div(class="ib","Only needed for replicates (multiple samples per condition). Skip for one sample per condition."),
      lapply(rv$samp_names,function(s)
        fluidRow(column(4,tags$b(s)," →"),
                 column(8,textInput(paste0("cond_",gsub("[^A-Za-z0-9]","_",s)),label=NULL,value=s)))
      ),
      actionButton("apply_conds","Apply Conditions",class="btn-sm btn-info"),br(),br()
    )
  })
  observeEvent(input$apply_conds,{
    req(rv$obj,rv$is_multi)
    tryCatch({
      sn<-rv$samp_names
      cmap<-setNames(sapply(sn,function(s){k<-paste0("cond_",gsub("[^A-Za-z0-9]","_",s));v<-input[[k]];if(is.null(v)||!nchar(trimws(v))) s else trimws(v)}),sn)
      if(!all(unique(rv$obj$sample)%in%names(cmap))){showNotification("Sample mismatch.",type="warning");return()}
      rv$obj$condition<-cmap[rv$obj$sample]
      if(!is.null(rv$obj_raw)) rv$obj_raw$condition<-cmap[rv$obj_raw$sample]
      if(!is.null(rv$obj_joined)) rv$obj_joined$condition<-cmap[rv$obj_joined$sample]
      rv$cond_map<-cmap
      showNotification(paste("Conditions:",paste(paste(names(cmap),cmap,sep="→"),collapse=", ")),type="message")
    },error=function(e) showNotification(paste("Condition error:",e$message),type="warning"))
  })

  # ---- QC ----
  observeEvent(input$apply_qc,{
    req(rv$obj_raw)
    withProgress(message="Applying QC...",value=0,{
      tryCatch({
        obj<-rv$obj_raw
        obj<-subset(obj,subset=nFeature_RNA>input$min_feat&nFeature_RNA<input$max_feat&percent.mt<input$max_mt)
        incProgress(0.3); obj<-NormalizeData(obj)
        incProgress(0.5); obj<-FindVariableFeatures(obj,nfeatures=2000)
        incProgress(0.7); obj<-ScaleData(obj)
        incProgress(0.85); obj<-RunPCA(obj,verbose=FALSE)
        rv$obj<-obj; rv$obj_filt<-obj; rv$qc_applied<-TRUE; rv$integ_done<-FALSE; rv$clust_done<-FALSE
        incProgress(1); showNotification(paste("QC done.",ncol(obj),"cells retained."),type="message")
      },error=function(e) showNotification(paste("QC error:",e$message),type="error"))
    })
  })

  # ---- INTEGRATION ----
  observeEvent(input$run_integ,{
    req(rv$obj_filt)
    meth<-input$integ_method
    if(!rv$is_multi&&meth!="none"){showNotification("Single sample — PCA only.",type="warning");meth<-"none"}
    withProgress(message="Running integration...",value=0,{
      tryCatch({
        obj<-rv$obj_filt
        if(meth=="sct_cca"){
          incProgress(0.1); obj<-SCTransform(obj,vars.to.regress="percent.mt",verbose=FALSE)
          incProgress(0.4); obj<-RunPCA(obj,verbose=FALSE)
          incProgress(0.6); set.seed(42); obj<-IntegrateLayers(obj,method=CCAIntegration,normalization.method="SCT",verbose=FALSE)
          ru<-"integrated.dr"
        } else if(meth=="cca"){
          incProgress(0.4); set.seed(42); obj<-IntegrateLayers(obj,method=CCAIntegration,orig.reduction="pca",new.reduction="integrated.cca",verbose=FALSE)
          ru<-"integrated.cca"
        } else if(meth=="harmony"){
          if(!requireNamespace("harmony",quietly=TRUE)) install.packages("harmony"); library(harmony)
          incProgress(0.4); set.seed(42); obj<-RunHarmony(obj,group.by.vars="sample",theta=input$harm_theta,plot_convergence=FALSE)
          ru<-"harmony"
        } else {ru<-"pca"}
        if(rv$is_multi){
          incProgress(0.75)
          set.seed(42); obj<-RunUMAP(obj,reduction="pca",dims=1:15,reduction.name="umap.unintegrated",seed.use=42,verbose=FALSE)
        }
        rv$obj<-obj; rv$integ_method<-meth; rv$reduc_use<-ru
        if(rv$is_multi && ru != "pca"){
          set.seed(42)
          obj <- RunUMAP(obj, reduction=ru, dims=1:input$n_dims, reduction.name="umap.integrated.preview", seed.use=42, verbose=FALSE)
          rv$obj <- obj
        }
        rv$harm_theta_used <- if(meth=="harmony") input$harm_theta else NA
        rv$integ_done<-TRUE; rv$clust_done<-FALSE; rv$markers_done<-FALSE
        incProgress(1); showNotification(paste("Integration done. Method:",meth),type="message")
      },error=function(e){showNotification(paste("Integration error:",e$message),type="error",duration=15);message(e$message)})
    })
  })

  # ---- CLUSTERING ----
  observeEvent(input$run_cluster,{
    req(rv$obj,rv$integ_done)
    withProgress(message="Clustering...",value=0,{
      tryCatch({
        obj<-rv$obj; dims<-1:input$n_dims; red<-rv$reduc_use
        incProgress(0.25); set.seed(42); obj<-FindNeighbors(obj,reduction=red,dims=dims,verbose=FALSE)
        incProgress(0.45); set.seed(42); obj<-FindClusters(obj,resolution=input$resolution,verbose=FALSE,random.seed=42)
        incProgress(0.65); set.seed(42); obj<-RunUMAP(obj,reduction=red,dims=dims,seed.use=42,verbose=FALSE)
        incProgress(0.85); obj[["RNA"]]<-JoinLayers(obj[["RNA"]])
        rv$obj<-obj; rv$clust_done<-TRUE
        rv$markers<-NULL; rv$markers_done<-FALSE; rv$cons_markers<-NULL; rv$cons_done<-FALSE
        rv$de_cache<-list(); rv$de_done<-FALSE; rv$manual_annot_done<-FALSE
        rv$sig_scores<-NULL; rv$sig_done<-FALSE; rv$cg_state<-NULL
        rv$pb_obj<-NULL; rv$pb_done<-FALSE; rv$gsea_res<-NULL; rv$gsea_done<-FALSE
        rv$ora_cache<-list(); rv$ora_res_up<-NULL; rv$ora_res_dn<-NULL; rv$ora_done<-FALSE
        rv$sig_gene_table<-NULL; rv$manual_annot_labels<-list()
        updateSelectInput(session,"comp_ident",choices=c("Cluster Number"="cluster"))
        incProgress(1); showNotification(paste("Clustering done.",length(levels(obj$seurat_clusters)),"clusters."),type="message")
      },error=function(e){showNotification(paste("Clustering error:",e$message),type="error",duration=15);message(e$message)})
    })
  })

  # ---- SUBCLUSTER: target cluster picker ----
  output$sub_target_clusters_ui <- renderUI({
    req(rv$obj,rv$clust_done)
    obj<-apply_ident(rv$obj,input$clust_ident%||%"cluster")
    ch<-nat_sort(unique(as.character(Idents(obj))))
    selectizeInput("sub_target_clusters","Cluster(s) to subcluster:",choices=ch,multiple=TRUE)
  })

  # ---- SUBCLUSTER: STAGE 1 — prep subset & PCA ----
  observeEvent(input$prep_subcluster,{
    req(rv$obj,input$sub_target_clusters)
    withProgress(message="Prepping subset...",value=0,{
      tryCatch({
        id_col<-get_id_col(input$clust_ident%||%"cluster")
        obj<-apply_ident(rv$obj,input$clust_ident%||%"cluster")
        sub<-subset(obj,idents=input$sub_target_clusters)
        incProgress(0.2); set.seed(42); sub<-FindVariableFeatures(sub,verbose=FALSE)
        incProgress(0.4); set.seed(42); sub<-ScaleData(sub,verbose=FALSE)
        incProgress(0.6); set.seed(42); sub<-RunPCA(sub,verbose=FALSE)
        if(rv$is_multi){
          incProgress(0.8)
          set.seed(42); sub<-RunUMAP(sub,reduction="pca",dims=1:min(15,ncol(sub[["pca"]])),reduction.name="umap.unintegrated",seed.use=42,verbose=FALSE)
        }
        rv$sub_obj<-sub; rv$sub_parent_clusters<-input$sub_target_clusters; rv$sub_parent_col<-id_col
        rv$subprep_done<-TRUE; rv$sub_integ_done<-FALSE; rv$sub_ru<-NULL; rv$sub_was_split<-FALSE; rv$subclust_done<-FALSE
        rv$sub_markers<-NULL; rv$sub_markers_done<-FALSE
        incProgress(1); showNotification(paste("Subset ready.",ncol(sub),"cells."),type="message")
      },error=function(e){showNotification(paste("Subcluster prep error:",e$message),type="error",duration=15);message(e$message)})
    })
  })

  # ---- SUBCLUSTER: stage-1 diagnostics (reuse global elbow/heatmap/pre-integration UMAP, pointed at rv$sub_obj) ----
  output$sub_diag_header <- renderUI({
    req(rv$subprep_done,rv$sub_parent_clusters)
    h4(paste0("Subset PCA Diagnostics — Cluster(s): ",paste(rv$sub_parent_clusters,collapse=", ")))
  })
  output$sub_elbow_ui <- renderUI({
    req(rv$sub_obj)
    ndims <- min(40,ncol(rv$sub_obj[["pca"]]))
    w_in <- max(8, ndims*0.15); h_in <- 5
    div(style="overflow:auto; max-width:100%; max-height:80vh;",
        plotOutput("sub_elbow", height=paste0(px_from_in(h_in),"px"), width=paste0(px_from_in(w_in),"px")))
  })
  output$sub_elbow    <- renderPlot({req(rv$sub_obj);ElbowPlot(rv$sub_obj,ndims=min(40,ncol(rv$sub_obj[["pca"]])))})
  output$sub_pca_heat_ui <- renderUI({
    req(rv$sub_obj)
    n <- min(input$n_heat%||%12,ncol(rv$sub_obj[["pca"]]))
    nrow_grid <- ceiling(n/3)
    w_in <- 12; h_in <- max(6, nrow_grid*3)
    div(style="overflow:auto; max-width:100%; max-height:80vh;",
        plotOutput("sub_pca_heat", height=paste0(px_from_in(h_in),"px"), width=paste0(px_from_in(w_in),"px")))
  })
  output$sub_pca_heat <- renderPlot({req(rv$sub_obj);DimHeatmap(rv$sub_obj,dims=1:min(input$n_heat%||%12,ncol(rv$sub_obj[["pca"]])),cells=500,balanced=TRUE,ncol=3)})
  sub_integ_view_reduction <- function(){
    if(!is.null(input$sub_integ_view) && input$sub_integ_view=="post" && !is.null(rv$sub_obj) && "umap.sub.integrated.preview"%in%names(rv$sub_obj@reductions)) "umap.sub.integrated.preview" else "umap.unintegrated"
  }
  sub_integ_view_title <- function() if(identical(sub_integ_view_reduction(),"umap.sub.integrated.preview")) "Subset — Post-integration (preview)" else "Subset — Unintegrated"
  output$sub_umap_pre_ui <- renderUI({
    req(rv$sub_obj,rv$is_multi)
    dd <- repel_dims(unique(as.character(rv$sub_obj$sample)))
    div(style="overflow:auto; max-width:100%; max-height:80vh;",
        plotOutput("sub_umap_pre", height=paste0(px_from_in(dd$h),"px"), width=paste0(px_from_in(dd$w),"px")))
  })
  output$sub_umap_pre <- renderPlot({
    req(rv$sub_obj,rv$is_multi,"umap.unintegrated"%in%names(rv$sub_obj@reductions))
    obj <- filter_by_samples(rv$sub_obj, input$active_samples)
    p <- umap_plot(obj,group.by="sample",pt.size=0.3,title=sub_integ_view_title(),reduction=sub_integ_view_reduction()) +
         fixed_ident_scale(rv$sub_obj$sample)
    grid::grid.draw(fix_panel_size(p))
  })

  # ---- SUBCLUSTER: STAGE 2 — integrate subset ----
  observeEvent(input$integ_subcluster,{
    req(rv$sub_obj,rv$subprep_done)
    withProgress(message="Integrating subset...",value=0,{
      tryCatch({
        sub<-rv$sub_obj
        meth<-switch(input$sub_integ_method,"Match global"=rv$integ_method,"None"="none",
                     "Harmony"="harmony","CCA"="cca","SCTransform+CCA"="sct_cca",input$sub_integ_method)
        if(!rv$is_multi&&meth!="none"){showNotification("Single sample — PCA only.",type="warning");meth<-"none"}
        theta_use<-if(identical(input$sub_integ_method,"Match global")) rv$harm_theta_used else input$sub_harm_theta
        if(identical(input$sub_integ_method,"Match global") && meth=="harmony" && is.na(theta_use)){
          showNotification("Original Harmony theta unknown (loaded from RDS) — using default theta of 2 for this subcluster run.",type="warning")
          theta_use<-2
        }
        was_split<-FALSE
        if(meth%in%c("cca","sct_cca")){
          lnames<-names(sub[["RNA"]]@layers)
          if(!any(grepl("\\.",lnames))){sub[["RNA"]]<-split(sub[["RNA"]],f=sub$sample); was_split<-TRUE}
        }
        incProgress(0.1)
        if(meth=="sct_cca"){
          incProgress(0.2); sub<-SCTransform(sub,vars.to.regress="percent.mt",verbose=FALSE)
          incProgress(0.4); sub<-RunPCA(sub,verbose=FALSE)
          incProgress(0.55); set.seed(42); sub<-IntegrateLayers(sub,method=CCAIntegration,normalization.method="SCT",verbose=FALSE)
          ru<-"integrated.dr"
        } else if(meth=="cca"){
          incProgress(0.4); set.seed(42); sub<-IntegrateLayers(sub,method=CCAIntegration,orig.reduction="pca",new.reduction="integrated.cca",verbose=FALSE)
          ru<-"integrated.cca"
        } else if(meth=="harmony"){
          if(!requireNamespace("harmony",quietly=TRUE)) install.packages("harmony"); library(harmony)
          incProgress(0.4); set.seed(42); sub<-RunHarmony(sub,group.by.vars="sample",theta=theta_use,plot_convergence=FALSE)
          ru<-"harmony"
        } else {ru<-"pca"}
        if(rv$is_multi && ru != "pca"){
          incProgress(0.8)
          set.seed(42)
          sub <- RunUMAP(sub, reduction=ru, dims=1:min(15,ncol(sub[[ru]])), reduction.name="umap.sub.integrated.preview", seed.use=42, verbose=FALSE)
        }
        rv$sub_obj<-sub; rv$sub_ru<-ru; rv$sub_was_split<-was_split; rv$sub_integ_done<-TRUE
        incProgress(1); showNotification(paste("Subset integration done. Method:",meth),type="message")
      },error=function(e){showNotification(paste("Subcluster integration error:",e$message),type="error",duration=15);message(e$message)})
    })
  })

  # ---- SUBCLUSTER: STAGE 3 — cluster ----
  observeEvent(input$run_subcluster,{
    req(rv$sub_obj,rv$sub_integ_done)
    withProgress(message="Running subclustering...",value=0,{
      tryCatch({
        sub<-rv$sub_obj; ru<-rv$sub_ru
        incProgress(0.3); set.seed(42); sub<-FindNeighbors(sub,reduction=ru,dims=1:input$sub_dims,verbose=FALSE)
        incProgress(0.55); set.seed(42); sub<-FindClusters(sub,resolution=input$sub_resolution,verbose=FALSE,random.seed=42)
        incProgress(0.8);  set.seed(42); sub<-RunUMAP(sub,reduction=ru,dims=1:input$sub_dims,seed.use=42,verbose=FALSE)
        if(rv$sub_was_split) sub[["RNA"]]<-JoinLayers(sub[["RNA"]])
        vals <- as.character(sub$seurat_clusters)
        lv <- nat_sort(unique(vals))
        prefix <- if(length(rv$sub_parent_clusters)==1) rv$sub_parent_clusters[1] else paste(sort(unique(rv$sub_parent_clusters)),collapse="+")
        sub$subcluster_preview_label <- paste0(prefix,"_",letters[match(vals,lv)])
        Idents(sub) <- sub$subcluster_preview_label
        rv$sub_obj<-sub; rv$subclust_done<-TRUE
        incProgress(1); showNotification(paste("Subclustering done.",length(levels(sub$seurat_clusters)),"subclusters."),type="message")
      },error=function(e){showNotification(paste("Subclustering error:",e$message),type="error",duration=15);message(e$message)})
    })
  })

  # ---- SUBCLUSTER: stage-2 result outputs ----
  output$sub_umap_ui <- renderUI({req(rv$sub_obj,rv$subclust_done)
    dd<-repel_dims(nat_sort(unique(as.character(rv$sub_obj$subcluster_preview_label))))
    plotOutput("sub_umap",height=paste0(px_from_in(dd$h),"px"),width=paste0(px_from_in(dd$w),"px"))})
  output$sub_umap <- renderPlot({
    req(rv$sub_obj,rv$subclust_done)
    obj <- filter_by_samples(rv$sub_obj, input$active_samples)
    p <- umap_plot(obj,label=TRUE,pt.size=0.4,title="Subcluster UMAP") +
         fixed_ident_scale(Idents(rv$sub_obj))
    grid::grid.draw(fix_panel_size(p))
  })
  output$sub_clust_dt <- renderDT({req(rv$sub_obj,rv$subclust_done)
    t<-table(rv$sub_obj$subcluster_preview_label)
    df<-data.frame(Subcluster=names(t),Count=as.numeric(t),Pct=round(as.numeric(t)/sum(t)*100,2))
    datatable(df,options=list(pageLength=30))})
  output$sub_clust_sample_dt <- renderDT({req(rv$sub_obj,rv$subclust_done)
    df <- data.frame(Subcluster=as.character(rv$sub_obj$subcluster_preview_label),
                     Sample=as.character(rv$sub_obj$sample))
    agg <- df %>% group_by(Subcluster,Sample) %>% summarise(Count=n(),.groups="drop") %>%
      group_by(Subcluster) %>% mutate(Pct=round(Count/sum(Count)*100,2)) %>% ungroup()
    datatable(agg,options=list(pageLength=30))})
  output$dl_sub_clust_dt <- downloadHandler(
    filename=function() paste0("subcluster_sizes_",Sys.Date(),".csv"),
    content=function(f){
      req(rv$sub_obj,rv$subclust_done)
      t<-table(rv$sub_obj$subcluster_preview_label)
      df<-data.frame(Subcluster=names(t),Count=as.numeric(t),Pct=round(as.numeric(t)/sum(t)*100,2))
      write.csv(df,f,row.names=FALSE)
    })
  output$dl_sub_clust_sample_dt <- downloadHandler(
    filename=function() paste0("subcluster_sizes_by_sample_",Sys.Date(),".csv"),
    content=function(f){
      req(rv$sub_obj,rv$subclust_done)
      df <- data.frame(Subcluster=as.character(rv$sub_obj$subcluster_preview_label),
                       Sample=as.character(rv$sub_obj$sample))
      agg <- df %>% group_by(Subcluster,Sample) %>% summarise(Count=n(),.groups="drop") %>%
        group_by(Subcluster) %>% mutate(Pct=round(Count/sum(Count)*100,2)) %>% ungroup()
      write.csv(agg,f,row.names=FALSE)
    })
  output$sub_clust_origin_dt <- renderDT({
    req(rv$sub_obj,rv$subclust_done)
    req(length(rv$sub_parent_clusters)>1)
    cells_sub<-colnames(rv$sub_obj)
    df <- data.frame(Subcluster=as.character(rv$sub_obj$subcluster_preview_label[cells_sub]),
                     Parent=as.character(rv$obj@meta.data[cells_sub,rv$sub_parent_col]))
    agg <- df %>% group_by(Subcluster,Parent) %>% summarise(Count=n(),.groups="drop") %>%
      group_by(Subcluster) %>% mutate(Pct=round(Count/sum(Count)*100,2)) %>% ungroup()
    datatable(agg,options=list(pageLength=30))})
  output$dl_sub_clust_origin_dt <- downloadHandler(
    filename=function() paste0("subcluster_sizes_by_parent_origin_",Sys.Date(),".csv"),
    content=function(f){
      req(rv$sub_obj,rv$subclust_done)
      req(length(rv$sub_parent_clusters)>1)
      cells_sub<-colnames(rv$sub_obj)
      df <- data.frame(Subcluster=as.character(rv$sub_obj$subcluster_preview_label[cells_sub]),
                       Parent=as.character(rv$obj@meta.data[cells_sub,rv$sub_parent_col]))
      agg <- df %>% group_by(Subcluster,Parent) %>% summarise(Count=n(),.groups="drop") %>%
        group_by(Subcluster) %>% mutate(Pct=round(Count/sum(Count)*100,2)) %>% ungroup()
      write.csv(agg,f,row.names=FALSE)
    })

  # ---- SUBCLUSTER: merge into main object ----
  observeEvent(input$merge_subcluster,{
    req(rv$obj,rv$sub_obj,rv$subclust_done)
    tryCatch({
      if(!"subcluster_labels"%in%colnames(rv$obj@meta.data))
        rv$obj@meta.data$subcluster_labels<-as.character(rv$obj@meta.data[[rv$sub_parent_col]])
      cells_sub<-colnames(rv$sub_obj)
      new_labels<-setNames(as.character(rv$sub_obj$subcluster_preview_label[cells_sub]),cells_sub)
      rv$obj@meta.data[cells_sub,"subcluster_labels"]<-unname(new_labels[cells_sub])
      rv$cg_state <- NULL
      rv$sig_scores <- NULL; rv$sig_done <- FALSE
      rv$markers <- NULL; rv$markers_done <- FALSE
      rv$cons_markers <- NULL; rv$cons_done <- FALSE
      rv$pb_obj <- NULL; rv$pb_done <- FALSE
      rv$de_cache <- list(); rv$de_done <- FALSE
      rv$gsea_res <- NULL; rv$gsea_done <- FALSE
      rv$ora_cache <- list(); rv$ora_res_up <- NULL; rv$ora_res_dn <- NULL; rv$ora_done <- FALSE
      showNotification(paste("Merged",length(unique(new_labels)),"new sub-label(s) across",length(rv$sub_parent_clusters),"parent cluster(s)."),type="message")
    },error=function(e){showNotification(paste("Merge error:",e$message),type="error",duration=15);message(e$message)})
  })

  # ---- SUBCLUSTER: discard ----
  observeEvent(input$discard_subcluster,{
    rv$sub_obj<-NULL; rv$subprep_done<-FALSE; rv$sub_integ_done<-FALSE; rv$sub_ru<-NULL; rv$sub_was_split<-FALSE; rv$subclust_done<-FALSE
    rv$sub_markers<-NULL; rv$sub_markers_done<-FALSE
  })

  # ---- FIND MARKERS ----
  observeEvent(input$find_markers,{
    req(rv$obj,rv$clust_done)
    withProgress(message="Finding markers...",value=0,{
      tryCatch({
        incProgress(0.3)
        obj <- apply_ident(rv$obj, input$clust_ident%||%"cluster")
        m<-FindAllMarkers(obj,only.pos=TRUE,min.pct=input$min_pct,logfc.threshold=input$min_lfc,verbose=FALSE)
        rv$markers<-m; rv$markers_done<-TRUE
        incProgress(1); showNotification("Markers done.",type="message")
      },error=function(e) showNotification(paste("Markers error:",e$message),type="error"))
    })
  })

  # ---- SUBCLUSTER: top markers (subclusters vs each other, scoped to rv$sub_obj) ----
  observeEvent(input$find_sub_markers,{
    req(rv$sub_obj,rv$subclust_done)
    withProgress(message="Finding subcluster markers...",value=0,{
      tryCatch({
        incProgress(0.3)
        obj <- rv$sub_obj
        Idents(obj) <- obj$subcluster_preview_label
        m<-FindAllMarkers(obj,only.pos=TRUE,min.pct=input$sub_min_pct,logfc.threshold=input$sub_min_lfc,verbose=FALSE)
        rv$sub_markers<-m; rv$sub_markers_done<-TRUE
        incProgress(1); showNotification("Subcluster markers done.",type="message")
      },error=function(e) showNotification(paste("Subcluster markers error:",e$message),type="error"))
    })
  })

  # ---- SCTYPE ----
  observeEvent(input$run_sctype,{
    req(rv$obj,rv$clust_done)
    withProgress(message="Running ScType...",value=0,{
      tryCatch({
        if(!requireNamespace("HGNChelper",quietly=TRUE)) install.packages("HGNChelper"); library(HGNChelper)
        source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
        source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
        gs<-gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx",input$tissue)
        sc<-as.matrix(rv$obj[["RNA"]]$scale.data)
        em<-sctype_score(scRNAseqData=sc,scaled=TRUE,gs=gs$gs_positive,gs2=gs$gs_negative)
        cl_res<-do.call("rbind",lapply(unique(rv$obj$seurat_clusters),function(cl){
          cells<-rownames(rv$obj@meta.data[rv$obj$seurat_clusters==cl,,drop=FALSE])
          es<-sort(rowSums(em[,cells,drop=FALSE]),decreasing=TRUE)
          head(data.frame(cluster=cl,type=names(es),scores=es,ncells=length(cells)),10)
        }))
        scs<-cl_res%>%group_by(cluster)%>%top_n(n=1,wt=scores)%>%ungroup()  # UNGROUP before mutations
        scs$type<-as.character(scs$type)
        scs$type[as.numeric(as.character(scs$scores))<scs$ncells/4]<-"Unknown"

        # Add marker_genes column — top genes that drove the ScType call
        scs$marker_genes <- sapply(seq_len(nrow(scs)), function(i) {
          ct <- scs$type[i]
          if(ct == "Unknown") return("—")
          pos <- gs$gs_positive[[ct]]
          if(is.null(pos) || !length(pos)) return("—")
          present <- pos[pos %in% rownames(sc)]
          if(!length(present)) return("—")
          scores_g <- rowMeans(sc[present, , drop=FALSE])
          paste(names(sort(scores_g, decreasing=TRUE))[seq_len(min(6, length(present)))], collapse=", ")
        })

        # Write back to seurat object
        rv$obj@meta.data$sctype_classification<-""
        for(j in unique(scs$cluster)){
          ct<-scs[scs$cluster==j,,drop=FALSE]
          rv$obj@meta.data$sctype_classification[rv$obj@meta.data$seurat_clusters==j]<-ct$type[1]
        }
        rv$sctype_res<-scs; rv$sctype_done<-TRUE

        ch<-c("Cluster Number"="cluster","ScType Label"="sctype")
        if("subcluster_labels"%in%colnames(rv$obj@meta.data)) ch<-c(ch,"Subcluster Label"="subcluster")
        updateSelectInput(session,"comp_ident",choices=ch)
        showNotification("ScType done.",type="message")
      },error=function(e) showNotification(paste("ScType error:",e$message),type="error"))
    })
  })

  # ---- SIGNATURE SCORING ----
  observeEvent(input$add_sig, {
    nm <- trimws(input$sig_name %||% "")
    genes_raw <- trimws(input$sig_genes %||% "")
    if(!nchar(nm)){showNotification("Enter a signature name.",type="warning");return()}
    if(!nchar(genes_raw)){showNotification("Enter at least one gene.",type="warning");return()}
    if(nm %in% names(rv$sig_sets)){showNotification(paste0("'",nm,"' already exists. Use a different name or Clear All."),type="warning");return()}
    genes <- trimws(unlist(strsplit(genes_raw,",")))
    genes <- genes[nchar(genes)>0]
    rv$sig_sets <- c(rv$sig_sets, setNames(list(genes), nm))
    updateTextInput(session,"sig_name",value="")
    updateTextInput(session,"sig_genes",value="")
    showNotification(paste("Added:",nm,"(",length(genes),"genes)"),type="message")
  })

  observeEvent(input$clear_sigs, {
    rv$sig_sets <- list(); rv$sig_done <- FALSE; rv$sig_scores <- NULL
  })

  output$sig_list_ui <- renderUI({
    if(!length(rv$sig_sets)) return(div(class="wb","No signatures added yet."))
    tagList(
      h5(paste(length(rv$sig_sets),"signature(s):")),
      lapply(names(rv$sig_sets), function(nm) {
        safe_nm <- gsub("[^A-Za-z0-9]","_",nm)
        fluidRow(
          column(3, tags$b(nm)),
          column(7, div(style="font-size:11px;color:#555", paste(rv$sig_sets[[nm]],collapse=", "))),
          column(2, actionButton(paste0("rm_sig_",safe_nm),"×",class="btn-danger btn-xs",
                     onclick=paste0("Shiny.setInputValue('rm_sig_key','",nm,"',{priority:'event'});")))
        )
      })
    )
  })
  observeEvent(input$load_sig_csv, {
    req(input$sig_csv, rv$obj)
    
    tryCatch({
      sig_df <- read.csv(input$sig_csv$datapath, stringsAsFactors = FALSE, check.names = FALSE)
      
      needed <- c("Signature", "Genes")
      if (!all(needed %in% colnames(sig_df))) {
        showNotification("CSV must contain columns named: Signature and Genes", type = "error", duration = 10)
        return()
      }
      
      sig_df <- sig_df[!is.na(sig_df$Signature) & !is.na(sig_df$Genes), , drop = FALSE]
      
      added <- 0
      skipped <- 0
      
      for (i in seq_len(nrow(sig_df))) {
        nm <- trimws(as.character(sig_df$Signature[i]))
        raw_genes <- as.character(sig_df$Genes[i])
        
        if (!nchar(nm) || !nchar(raw_genes)) {
          skipped <- skipped + 1
          next
        }
        
        genes <- resolve_genes(raw_genes, rv$obj)
        
        if (length(genes) > 0) {
          rv$sig_sets[[nm]] <- genes
          added <- added + 1
        } else {
          skipped <- skipped + 1
        }
      }
      
      rv$sig_done <- FALSE
      rv$sig_scores <- NULL
      rv$sig_gene_table <- NULL
      
      showNotification(
        paste0("Loaded ", added, " signature(s). Skipped ", skipped, "."),
        type = "message",
        duration = 8
      )
      
    }, error = function(e) {
      showNotification(paste("Signature CSV error:", e$message), type = "error", duration = 10)
    })
  })
  observeEvent(input$rm_sig_key, {
    key <- input$rm_sig_key
    rv$sig_sets[[key]] <- NULL
  })

  observeEvent(input$run_sig, {
    req(rv$obj, rv$clust_done)
    if(!length(rv$sig_sets)){showNotification("Add at least one signature first.",type="warning");return()}
    withProgress(message="Scoring signatures...",value=0,{
      tryCatch({
        obj <- rv$obj; all_genes <- rownames(obj)
        id_col <- get_id_col(input$clust_ident%||%"cluster")
        scores_list <- list(); missing_msgs <- character(0)
        for(nm in names(rv$sig_sets)){
          incProgress(1/length(rv$sig_sets),detail=nm)
          raw_g <- rv$sig_sets[[nm]]
          valid  <- raw_g[raw_g %in% all_genes]
          miss   <- setdiff(raw_g, all_genes)
          if(length(miss)) missing_msgs <- c(missing_msgs, paste0(nm,":[",paste(miss,collapse=","),"]"))
          if(length(valid)<2){
            showNotification(paste0("'",nm,"' has <2 valid genes — skipped."),type="warning"); next
          }
          ctrl_n <- min(100, length(valid))
          safe_nm <- gsub("[^A-Za-z0-9]","_",nm)
          col_nm  <- paste0("sig_",safe_nm,"1")
          obj <- AddModuleScore(obj,features=list(valid),name=paste0("sig_",safe_nm),ctrl=ctrl_n,assay="RNA")
          scores_list[[nm]] <- tapply(obj@meta.data[[col_nm]], obj@meta.data[[id_col]], mean)
        }
        if(length(missing_msgs))
          showNotification(paste("Genes not found (skipped):",paste(missing_msgs,collapse="; ")),type="warning",duration=12)
        if(!length(scores_list)){showNotification("No signatures scored.",type="warning");return()}
        cls_ord <- nat_sort(unique(as.character(obj@meta.data[[id_col]])))
        mat <- do.call(rbind, lapply(scores_list, function(x) as.numeric(x[cls_ord])))
        colnames(mat) <- cls_ord; rownames(mat) <- names(scores_list)
        rv$sig_scores <- mat; rv$sig_done <- TRUE
        # Build long-format gene inspector table
        tryCatch({
          expr_mat <- GetAssayData(obj, assay="RNA", layer="data")
          gene_rows <- list()
          for(sig_nm in names(rv$sig_sets)) {
            raw_g <- rv$sig_sets[[sig_nm]]
            found_g <- raw_g[raw_g %in% rownames(expr_mat)]
            not_found <- setdiff(raw_g, rownames(expr_mat))
            for(g in found_g) {
              gene_rows[[length(gene_rows)+1]] <- data.frame(
                Signature=sig_nm, Gene=g, Found=TRUE,
                Mean_Expression=round(mean(expr_mat[g,],na.rm=TRUE),4), stringsAsFactors=FALSE)
            }
            for(g in not_found) {
              gene_rows[[length(gene_rows)+1]] <- data.frame(
                Signature=sig_nm, Gene=g, Found=FALSE,
                Mean_Expression=NA_real_, stringsAsFactors=FALSE)
            }
          }
          if(length(gene_rows)) {
            gt <- do.call(rbind, gene_rows)
            found_vals <- gt$Mean_Expression[gt$Found & !is.na(gt$Mean_Expression)]
            if(length(found_vals) >= 3) {
              qt <- quantile(found_vals, c(1/3, 2/3), na.rm=TRUE)
              gt$Contribution <- ifelse(!gt$Found, "Absent",
                                 ifelse(is.na(gt$Mean_Expression), "Absent",
                                 ifelse(gt$Mean_Expression >= qt[2], "High",
                                 ifelse(gt$Mean_Expression >= qt[1], "Medium", "Low"))))
            } else {
              gt$Contribution <- ifelse(!gt$Found, "Absent", "High")
            }
            rv$sig_gene_table <- gt
          }
        }, error=function(e) message("sig_gene_table build error: ", e$message))
        showNotification(paste("Done:",nrow(mat),"signatures across",ncol(mat),"clusters."),type="message")
      },error=function(e) showNotification(paste("Scoring error:",e$message),type="error"))
    })
  })

  output$sig_heat_ui <- renderUI({
    req(rv$sig_done, rv$sig_scores)
    mat <- rv$sig_scores
    w_in <- max(10, ncol(mat)*0.45 + max_label_len(colnames(mat))*0.12)
    h_in <- max(8, nrow(mat)*0.45 + max_label_len(rownames(mat))*0.12)
    plotOutput("sig_heat", height=paste0(px_from_in(h_in),"px"), width=paste0(px_from_in(w_in),"px"))
  })
  output$sig_heat <- renderPlot({
    req(rv$sig_done, rv$sig_scores)
    mat <- rv$sig_scores
    zmat <- if(nrow(mat)>1) t(scale(t(mat))) else mat
    zmat[!is.finite(zmat)] <- 0
    col_margin <- max(8, max_label_len(colnames(zmat))*0.6)
    row_margin <- max(10, max_label_len(rownames(zmat))*0.6)
    heatmap(zmat, scale="none",
            col=colorRampPalette(c("#3498db","white","#e74c3c"))(100),
            margins=c(col_margin, row_margin), main="Signature Scores — z-scored per signature",
            cexRow=1.0, cexCol=1.0)
  }, res=72)

  output$sig_tbl <- renderDT({
    req(rv$sig_done, rv$sig_scores)
    df <- cbind(Signature=rownames(rv$sig_scores), as.data.frame(round(rv$sig_scores,4)))
    datatable(df, rownames=FALSE, options=list(pageLength=20, scrollX=TRUE))
  })

  output$sig_view_sel_ui <- renderUI({
    req(rv$sig_done, rv$sig_gene_table)
    selectInput("sig_view_select","View signature:",choices=c("All",names(rv$sig_sets)),selected="All")
  })
  output$sig_gene_dt <- renderDT({
    req(rv$sig_done, rv$sig_gene_table)
    df <- rv$sig_gene_table
    if(!is.null(input$sig_view_select) && input$sig_view_select != "All")
      df <- df[df$Signature == input$sig_view_select, , drop=FALSE]
    datatable(df, rownames=FALSE, filter="top", options=list(pageLength=20, scrollX=TRUE))
  })
  output$dl_sig_gene_csv <- downloadHandler(
    filename=function() paste0("sig_gene_contributions_",Sys.Date(),".csv"),
    content=function(f){ req(rv$sig_gene_table); write.csv(rv$sig_gene_table, f, row.names=FALSE) }
  )

  # ---- MANUAL ANNOTATION ----
  mannot_refresh <- reactiveVal(0)
  observeEvent(input$refresh_mannot, mannot_refresh(mannot_refresh()+1))

  output$manual_annot_ui <- renderUI({
    mannot_refresh()
    req(rv$clust_done)

    id_col <- get_id_col(input$clust_ident%||%"cluster")
    cls <- nat_sort(unique(as.character(rv$obj@meta.data[[id_col]])))

    ctm <- if(isTRUE(rv$sctype_done) && !is.null(rv$sctype_res))
      setNames(as.character(rv$sctype_res$type), as.character(rv$sctype_res$cluster))
    else setNames(rep("", length(cls)), cls)
    
    mkr_top <- if(isTRUE(rv$markers_done) && !is.null(rv$markers)) {
      rv$markers %>%
        group_by(cluster) %>%
        top_n(n=5, wt=avg_log2FC) %>%
        summarise(top5=paste(gene, collapse=", "), .groups="drop")
    } else data.frame(cluster=cls, top5="")
    
    tagList(lapply(cls, function(cl) {
      
      sctype_hint <- ctm[[cl]] %||% ""
      
      top5 <- mkr_top$top5[match(cl, as.character(mkr_top$cluster))]
      top5 <- top5 %||% ""
      
      sig_hint <- if(isTRUE(rv$sig_done) &&
                     !is.null(rv$sig_scores) &&
                     cl %in% colnames(rv$sig_scores)) {
        
        vals <- rv$sig_scores[, cl]
        vals <- sort(vals, decreasing = TRUE)
        paste(names(vals)[seq_len(min(3, length(vals)))], collapse = ", ")
        
      } else {
        "run Signature Scoring first"
      }
      
      wellPanel(
        fluidRow(
          column(
            3,
            tags$b(paste("Cluster", cl)),
            div(style="font-size:11px;color:#777",
                paste0("ScType: ", sctype_hint)),
            div(style="font-size:10px;color:#999",
                paste0("Top genes: ", top5)),
            div(style="font-size:10px;color:#27ae60",
                paste0("Top sigs: ", sig_hint))
          ),
          column(
            9,
            textInput(
              paste0("mannot_", cl),
              "Manual label:",
              value = {
                parts <- strsplit(cl, "_", fixed=TRUE)[[1]]
                ancestor_hint <- NULL
                if(length(parts) > 1){
                  for(i in (length(parts)-1):1){
                    anc <- paste(parts[1:i], collapse="_")
                    if(!is.null(rv$manual_annot_labels[[anc]])){
                      suffix <- paste(parts[(i+1):length(parts)], collapse="_")
                      ancestor_hint <- paste0(rv$manual_annot_labels[[anc]], "_", suffix)
                      break
                    }
                  }
                }
                # subcluster_labels values are always already-meaningful text -- either
                # inherited verbatim from whatever identity was active before subclustering,
                # or a generated "<parent>_<letter>" combo -- never a bare numeric id, so
                # use `cl` as-is instead of prepending "Cluster_".
                if(is.null(ancestor_hint) && id_col=="subcluster_labels") ancestor_hint <- cl
                rv$manual_annot_labels[[cl]] %||% ancestor_hint %||% sctype_hint %||% paste0("Cluster_", cl)
              },
              placeholder = "Type annotation..."
            )
          )
        )
      )
    }))
  })

  observeEvent(input$apply_manual_annot, {
    req(rv$obj, rv$clust_done)
    tryCatch({
      id_col <- get_id_col(input$clust_ident%||%"cluster")
      cls <- nat_sort(unique(as.character(rv$obj@meta.data[[id_col]])))
      lm <- sapply(cls, function(cl) {
        v <- input[[paste0("mannot_", cl)]]
        if(is.null(v) || !nchar(trimws(v))) paste0("Cluster_", cl) else trimws(v)
      })
      names(lm) <- cls
      rv$obj@meta.data$manual_annotation <- lm[as.character(rv$obj@meta.data[[id_col]])]
      rv$manual_annot_labels <- as.list(lm)
      rv$manual_annot_done <- TRUE
      updateSelectInput(session, "comp_ident", choices=get_ident_ch())
      showNotification(paste("Manual labels applied.", length(unique(rv$obj$manual_annotation)), "unique labels."), type="message")
    }, error=function(e) showNotification(paste("Manual annotation error:", e$message), type="error"))
  })

  observeEvent(input$reset_manual_annot, {
    req(rv$clust_done)
    id_col <- get_id_col(input$clust_ident%||%"cluster")
    cls <- nat_sort(unique(as.character(rv$obj@meta.data[[id_col]])))
    ctm <- if(isTRUE(rv$sctype_done) && !is.null(rv$sctype_res))
      setNames(as.character(rv$sctype_res$type), as.character(rv$sctype_res$cluster))
    else setNames(rep("", length(cls)), cls)
    for(cl in cls)
      updateTextInput(session, paste0("mannot_", cl), value=ctm[cl] %||% "")
    showNotification("Text inputs reset to ScType suggestions.", type="message")
  })

  output$manual_annot_status <- renderPrint({
    req(rv$manual_annot_done, rv$obj)
    cat("Manual annotation applied.\n")
    t <- table(rv$obj$manual_annotation)
    for(n in names(t)) cat(" ", n, ":", t[[n]], "cells\n")
  })

  # ---- COMPARATIVE ----
  observeEvent(input$run_cons,{
    req(rv$obj,rv$clust_done,rv$is_multi,input$cons_clust)
    withProgress(message="Conserved markers...",value=0,{
      tryCatch({
        obj<-rv$obj; id_col<-get_id_col(input$comp_ident); Idents(obj)<-obj@meta.data[[id_col]]; cl<-as.character(input$cons_clust)
        incProgress(0.3)
        cm<-FindConservedMarkers(obj,ident.1=cl,grouping.var="sample",assay="RNA",verbose=FALSE)
        rv$cons_markers<-cbind(gene=rownames(cm),cm); rv$cons_done<-TRUE
        incProgress(1); showNotification("Conserved markers done.",type="message")
      },error=function(e){showNotification(paste("Conserved markers error:",e$message),type="error");message(e$message)})
    })
  })

  observeEvent(input$run_pb,{
    req(rv$obj,rv$clust_done,rv$is_multi)
    withProgress(message="Pseudobulk...",value=0,{
      tryCatch({
        id_col <- get_id_col(input$comp_ident)
        incProgress(0.4)
        # Seurat v5 AggregateExpression sanitises group names via make.names()
        # We store a lookup from original label -> sanitised colname
        agg <- AggregateExpression(rv$obj, group.by=c(id_col,"sample"),
                                   return.seurat=TRUE)
        rv$pb_obj      <- agg
        rv$pb_done     <- TRUE
        rv$pb_id_col   <- id_col
        incProgress(1)
        showNotification(
          paste0("Pseudobulk done. ", length(colnames(agg)), " columns created."),
          type="message")
      },error=function(e) showNotification(paste("Pseudobulk error:",e$message),type="error"))
    })
  })

  observeEvent(input$run_de,{
    req(rv$obj,rv$clust_done,rv$is_multi,input$de_clust)
    s1<-get_s1(); s2<-get_s2()
    withProgress(message="Running DE...",value=0,{
      tryCatch({
        obj<-rv$obj; id_col<-get_id_col(input$comp_ident)
        obj$cs<-paste(obj@meta.data[[id_col]],obj$sample,sep="_"); Idents(obj)<-"cs"
        i1<-paste(as.character(input$de_clust),s1,sep="_")
        i2<-paste(as.character(input$de_clust),s2,sep="_")
        ids<-levels(Idents(obj))
        if(!i1%in%ids||!i2%in%ids){
          showNotification(paste("Identity not found.\nAvailable:",paste(ids,collapse=", ")),type="warning"); return()
        }
        incProgress(0.4); de<-FindMarkers(obj,ident.1=i1,ident.2=i2,verbose=FALSE)
        de$gene<-rownames(de); de<-de[,c("gene",setdiff(colnames(de),"gene"))]
        rv$de_res<-de; rv$de_done<-TRUE; rv$de_s1<-s1; rv$de_s2<-s2
        # cache alongside the single-result slots above (not a replacement) —
        # foundational plumbing for upcoming GSEA/ORA tabs that need to look
        # back at more than just the most recent DE run
        de_key<-paste(id_col,as.character(input$de_clust),s1,s2,sep="|")
        rv$de_cache[[de_key]]<-list(de_res=de,s1=s1,s2=s2,
                                     cluster=as.character(input$de_clust),
                                     id_col=id_col,timestamp=Sys.time())
        incProgress(1); showNotification("DE done.",type="message")
      },error=function(e) showNotification(paste("DE error:",e$message),type="error"))
    })
  })

  # ---- GSEA ----
  # species+collection -> pathways_list/term2gene. Plain reactive (lazy) —
  # msigdbr() is only ever actually called when something reads this, which
  # only happens inside the run_gsea observer below, so changing the species/
  # collection dropdowns alone never triggers a fetch. Reused as-is by a
  # future ORA tab.
  gsea_pathways <- reactive({
    req(input$gsea_species, input$gsea_collection)
    fetch_msig(input$gsea_species, input$gsea_collection)
  })

  ora_pathways <- reactive({
    req(input$ora_species, input$ora_collection)
    fetch_msig(input$ora_species, input$ora_collection)
  })

  output$has_de_cache <- reactive(length(rv$de_cache)>0)
  outputOptions(output,"has_de_cache",suspendWhenHidden=FALSE)

  output$gsea_comparison_ui <- renderUI({
    req(length(rv$de_cache)>0)
    labs<-sapply(rv$de_cache,function(x) paste0(x$cluster,": ",x$s1," vs ",x$s2))
    selectInput("gsea_comparison","DE comparison (from cache):",choices=setNames(names(rv$de_cache),labs))
  })

  observeEvent(input$run_gsea,{
    req(input$gsea_comparison,input$gsea_species,input$gsea_collection)
    ck<-rv$de_cache[[input$gsea_comparison]]
    if(is.null(ck)){showNotification("Selected comparison not found in cache — re-run Condition DE.",type="warning");return()}
    withProgress(message="Running GSEA...",value=0,{
      tryCatch({
        incProgress(0.2,detail="Fetching pathways...")
        pl<-gsea_pathways()$pathways_list
        if(!length(pl)){showNotification("No pathways found for this species/collection.",type="warning");return()}
        incProgress(0.5,detail="Ranking genes...")
        stats<-build_gsea_stats(ck$de_res)
        if(length(stats)<2){showNotification("Not enough finite-ranked genes for GSEA.",type="warning");return()}
        incProgress(0.7,detail="Running fgsea...")
        fres<-tryCatch(fgseaMultilevel(pathways=pl,stats=stats),error=function(e) NULL)
        if(is.null(fres)||nrow(fres)==0){showNotification("fgsea returned no results.",type="warning");return()}
        fres<-fres%>%arrange(padj)%>%
          mutate(direction=ifelse(NES>0,paste("Higher in",ck$s1),paste("Higher in",ck$s2)))
        rv$gsea_res<-fres; rv$gsea_done<-TRUE; rv$gsea_s1<-ck$s1; rv$gsea_s2<-ck$s2
        rv$gsea_label<-paste0(ck$cluster,": ",ck$s1," vs ",ck$s2)
        rv$gsea_pl<-pl; rv$gsea_stats<-stats
        incProgress(1); showNotification("GSEA done.",type="message")
      },error=function(e) showNotification(paste("GSEA error:",e$message),type="error",duration=15))
    })
  })

  # ---- ORA ----
  output$ora_comparison_ui <- renderUI({
    req(length(rv$de_cache)>0)
    labs <- sapply(rv$de_cache, function(x) paste0(x$cluster,": ",x$s1," vs ",x$s2))
    selectInput("ora_comparison","DE comparison (from cache):",choices=setNames(names(rv$de_cache),labs))
  })

  filter_top_ora <- function(res_df, p_cut, topn) {
    if(is.null(res_df)||!nrow(res_df)) return(NULL)
    sig <- res_df[is.finite(res_df$p.adjust) & res_df$p.adjust < p_cut,,drop=FALSE]
    if(!nrow(sig)) return(NULL)
    sig[order(sig$p.adjust),,drop=FALSE][seq_len(min(topn,nrow(sig))),,drop=FALSE]
  }

  observeEvent(input$run_ora,{
    req(input$ora_comparison, input$ora_species, input$ora_collection)
    ck <- rv$de_cache[[input$ora_comparison]]
    if(is.null(ck)){showNotification("Selected comparison not found in cache.",type="warning");return()}
    withProgress(message="Running ORA...",value=0,{
      tryCatch({
        incProgress(0.2,detail="Fetching pathways...")
        t2g <- ora_pathways()$term2gene
        if(is.null(t2g)||nrow(t2g)==0){showNotification("No pathways found for this species/collection.",type="warning");return()}
        incProgress(0.4,detail="Splitting genes...")
        de       <- ck$de_res
        p_cut    <- as.numeric(input$ora_pcut %||% 0.05)
        lfc_cut  <- as.numeric(input$ora_lfc  %||% 0.5)
        topn     <- as.integer(input$ora_topn %||% 20)
        finite_ok <- is.finite(de$avg_log2FC) & is.finite(de$p_val_adj)
        sig_up   <- de$gene[finite_ok & de$avg_log2FC >  lfc_cut & de$p_val_adj < p_cut]
        sig_dn   <- de$gene[finite_ok & de$avg_log2FC < -lfc_cut & de$p_val_adj < p_cut]
        incProgress(0.6,detail="Running enricher()...")
        raw_up <- if(length(sig_up)>=2) tryCatch(as.data.frame(enricher(sig_up,TERM2GENE=t2g)),error=function(e) NULL) else NULL
        raw_dn <- if(length(sig_dn)>=2) tryCatch(as.data.frame(enricher(sig_dn,TERM2GENE=t2g)),error=function(e) NULL) else NULL
        rv$ora_res_up <- filter_top_ora(raw_up, p_cut, topn)
        rv$ora_res_dn <- filter_top_ora(raw_dn, p_cut, topn)
        rv$ora_done   <- TRUE
        rv$ora_s1     <- ck$s1; rv$ora_s2 <- ck$s2
        rv$ora_label  <- paste0(ck$cluster,": ",ck$s1," vs ",ck$s2)
        ora_key <- paste(ck$id_col,ck$cluster,ck$s1,ck$s2,sep="|")
        rv$ora_cache[[ora_key]] <- list(res_up=rv$ora_res_up,res_dn=rv$ora_res_dn,
                                         s1=ck$s1,s2=ck$s2,label=rv$ora_label)
        incProgress(1)
        n_up <- if(!is.null(rv$ora_res_up)) nrow(rv$ora_res_up) else 0
        n_dn <- if(!is.null(rv$ora_res_dn)) nrow(rv$ora_res_dn) else 0
        showNotification(paste0("ORA done. ",n_up," terms higher in ",ck$s1,
                                ", ",n_dn," terms higher in ",ck$s2,"."),type="message")
      },error=function(e) showNotification(paste("ORA error:",e$message),type="error",duration=15))
    })
  })

  observeEvent(input$plot_comp,{
    req(rv$obj,rv$clust_done,rv$is_multi,input$comp_genes)
    gns<-resolve_genes(input$comp_genes,rv$obj)
    if(!length(gns)){showNotification("No genes found.",type="warning");return()}
    obj<-apply_ident(rv$obj,input$comp_ident%||%"cluster")
    id_col<-get_id_col(input$comp_ident); ns<-length(rv$samp_names); ng<-length(gns)
    feat_ht <- max(600, ng*280)
    output$cg_feat_ui <- renderUI(plotOutput("cg_feat_inner", height=paste0(feat_ht,"px")))
    output$cg_feat_inner <- renderPlot({
      lims <- umap_full_xylim(rv$obj)
      obj_feat <- filter_by_samples(obj, input$active_samples)
      ns_feat <- length(unique(as.character(obj_feat$sample)))
      FeaturePlot(obj_feat,features=gns,split.by="sample",reduction="umap",max.cutoff=3,pt.size=0.2,ncol=ns_feat) &
        coord_cartesian(xlim=lims$x, ylim=lims$y)
    })
    # cached so the global/per-cluster scope toggle (cg_vln_scope/cg_dot_scope)
    # can re-render + re-size reactively without requiring another Plot click
    rv$cg_state <- list(obj=obj, gns=gns, id_col=id_col, ns=ns)
  })

  output$cg_vln_clust_ui <- renderUI({
    req(rv$cg_state)
    ch<-nat_sort(unique(as.character(rv$cg_state$obj@meta.data[[rv$cg_state$id_col]])))
    selectInput("cg_vln_clust","Cluster:",ch,selected=ch[1])
  })
  output$cg_dot_clust_ui <- renderUI({
    req(rv$cg_state)
    ch<-nat_sort(unique(as.character(rv$cg_state$obj@meta.data[[rv$cg_state$id_col]])))
    selectInput("cg_dot_clust","Cluster:",ch,selected=ch[1])
  })
  output$cg_vln_clusters_sel_ui <- renderUI({
    req(rv$cg_state)
    ch<-nat_sort(unique(as.character(rv$cg_state$obj@meta.data[[rv$cg_state$id_col]])))
    selectizeInput("cg_vln_clusters_sel","Clusters:",choices=ch,multiple=TRUE)
  })
  output$cg_dot_clusters_sel_ui <- renderUI({
    req(rv$cg_state)
    ch<-nat_sort(unique(as.character(rv$cg_state$obj@meta.data[[rv$cg_state$id_col]])))
    selectizeInput("cg_dot_clusters_sel","Clusters:",choices=ch,multiple=TRUE)
  })

  output$cg_vln_ui <- renderUI({
    req(rv$cg_state)
    st <- rv$cg_state
    
    clust_arg <- if(
      identical(input$cg_vln_scope %||% "global", "selected")
    ) {
      input$cg_vln_clusters_sel
    } else {
      input$cg_vln_clust
    }
    
    obj_use <- cg_scope_obj(
      st$obj,
      st$id_col,
      input$cg_vln_scope %||% "global",
      clust_arg
    )
    
    id_labels <- unique(
      as.character(obj_use@meta.data[[st$id_col]])
    )
    
    # Match the more spacious architecture used by the independent
    # Gene Expression violin plot: roughly 200 px per stacked gene.
    # When scope isn't "percluster", ggplot's own rotated (45deg) identity
    # tick labels consume their own separate auto-sized space in addition
    # to cg_vln_bottom_margin's explicit plot.margin allocation -- without
    # budgeting for that here too, long identity labels can starve the
    # panel down to near-zero height (see the matching comment on h_final
    # in dl_cg_vln for the full diagnosis). "percluster" scope resets the
    # label angle to horizontal via cg_vln_apply_single, so it never needs
    # this extra budget.
    ht <- max(500, length(st$gns) * 200) +
      cg_vln_bottom_margin(id_labels) +
      if (!identical(input$cg_vln_scope %||% "global", "percluster")) {
        cg_vln_bottom_margin(id_labels)
      } else {
        0
      }

    is_single <- identical(
      input$cg_vln_scope %||% "global",
      "percluster"
    )
    
    if (is_single) {
      div(
        style = "overflow:auto; max-width:100%; max-height:80vh;",
        plotOutput(
          "cg_vln_inner",
          height = paste0(ht, "px"),
          width = paste0(
            px_from_in(cg_vln_single_width_in),
            "px"
          )
        )
      )
    } else {
      w_use <- cg_vln_width_in(id_labels)
      # Floor in px assuming Shiny's default 72dpi renderPlot device (1pt==1px
      # here), matching cg_vln_min_canvas_in's inch-based calibration -- see
      # cg_vln_width_in's comment for why this floor exists.
      ht <- max(ht, cg_vln_min_canvas_in(id_labels)$h * 72)

      div(
        style = "overflow:auto; max-width:100%; max-height:80vh;",
        plotOutput(
          "cg_vln_inner",
          height = paste0(ht, "px"),
          width = paste0(px_from_in(w_use), "px")
        )
      )
    }
  })
  output$cg_vln_inner <- renderPlot({
    req(rv$cg_state)
    st <- rv$cg_state
    
    clust_arg <- if(
      identical(input$cg_vln_scope %||% "global", "selected")
    ) {
      input$cg_vln_clusters_sel
    } else {
      input$cg_vln_clust
    }
    
    obj_use <- cg_scope_obj(
      st$obj,
      st$id_col,
      input$cg_vln_scope %||% "global",
      clust_arg
    )
    obj_use <- filter_by_samples(obj_use, input$active_samples)

    id_labels <- unique(
      as.character(obj_use@meta.data[[st$id_col]])
    )

    gns_use <- nonconstant_genes(obj_use, st$gns)
    
    is_single <- identical(
      input$cg_vln_scope %||% "global",
      "percluster"
    )
    
    if (!length(gns_use)) {
      return(
        ggplot() +
          theme_void() +
          annotate(
            "text",
            x = 0,
            y = 0,
            label = paste(
              "All selected genes are constant",
              "(no variation) in this cluster."
            )
          )
      )
    }
    
    side_margin <- cg_vln_side_margin(id_labels)
    
    if (length(gns_use) == 1) {
      # Seurat's stack=TRUE path requires at least two features.
      p <- VlnPlot(
        obj_use,
        features = gns_use,
        group.by = st$id_col,
        split.by = "sample",
        pt.size = 0
      ) +
        geom_point(
          position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.9),
          size = 0.3,
          alpha = 0.6
        ) &
        theme(
          plot.margin = margin(
            t = 60,
            r = side_margin,
            b = cg_vln_bottom_margin(id_labels),
            l = side_margin,
            unit = "pt"
          )
        )
    } else {
      p <- VlnPlot(
        obj_use,
        features = gns_use,
        group.by = st$id_col,
        split.by = "sample",
        stack = TRUE,
        flip = TRUE,
        pt.size = 0
      ) +
        geom_point(
          position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.9),
          size = 0.05,
          alpha = 0.15
        ) &
        theme(
          plot.margin = margin(
            t = 60,
            r = max(
              side_margin,
              max_label_len(gns_use) * 4
            ),
            b = cg_vln_bottom_margin(id_labels),
            l = side_margin,
            unit = "pt"
          ),
          panel.spacing = unit(0.6, "cm")
        )
    }
    
    # Shared styling for both one-gene and stacked violin plots.
    # Keeping the y expansion here prevents either branch from using
    # Seurat's tight default scale.
    p <- p &
      min_text_theme &
      scale_y_continuous(
        expand = expansion(mult = c(0.05, 0.15))
      ) &
      theme(
        axis.text.x = element_text(
          size = if (length(id_labels) <= 28) {
            13
          } else {
            pmax(
              9,
              13 - (length(id_labels) - 28) * 0.3
            )
          }
        ),
        axis.text.y = element_text(
          size = if (length(id_labels) <= 28) {
            13
          } else {
            pmax(
              9,
              13 - (length(id_labels) - 28) * 0.3
            )
          }
        ),
        axis.title = element_text(size = 12)
      )
    
    cg_vln_apply_single(
      p,
      is_single,
      id_labels[1]
    )
  })

  output$cg_dot_ui <- renderUI({
    req(rv$cg_state); st<-rv$cg_state
    clust_arg<-if(identical(input$cg_dot_scope%||%"global","selected")) input$cg_dot_clusters_sel else input$cg_dot_clust
    obj_use<-cg_scope_obj(st$obj,st$id_col,input$cg_dot_scope%||%"global",clust_arg)
    id_labels<-unique(as.character(obj_use@meta.data[[st$id_col]]))
    dd <- dotplot_dims(length(st$gns), length(id_labels)*st$ns, id_labels)
    plotOutput("cg_dot", height=paste0(px_from_in(dd$h),"px"),
                          width=paste0(px_from_in(dd$w),"px"))
  })
  output$cg_dot<-renderPlot({
    req(rv$cg_state); st<-rv$cg_state
    clust_arg<-if(identical(input$cg_dot_scope%||%"global","selected")) input$cg_dot_clusters_sel else input$cg_dot_clust
    obj_use<-cg_scope_obj(st$obj,st$id_col,input$cg_dot_scope%||%"global",clust_arg)
    obj_use<-filter_by_samples(obj_use,input$active_samples)
    id_labels<-unique(as.character(obj_use@meta.data[[st$id_col]]))
    samp_names_use<-sort(unique(as.character(obj_use$sample)))
    ns_use<-length(samp_names_use)
    cols_u<-rep(c("blue","red","green","orange","purple"),ns_use)[1:ns_use]
    p<-DotPlot(obj_use,features=st$gns,split.by="sample",cols=cols_u,dot.scale=8)+RotatedAxis()&
      theme(plot.margin=margin(t=10,r=10,b=max(60,max_label_len(id_labels)*3),l=10,unit="pt"))&min_text_theme&
      theme(axis.text.x=element_text(size=if(length(id_labels)<=28) 13 else pmax(9, 13-(length(id_labels)-28)*0.3)), axis.text.y=element_text(size=if(length(id_labels)<=28) 13 else pmax(9, 13-(length(id_labels)-28)*0.3)), axis.title=element_text(size=12))
    with_dot_legend(p,samp_names_use,cols_u)
  })

  observeEvent(input$plot_gene,{
    req(rv$obj,rv$clust_done,input$gene_in)
    gns<-resolve_genes(input$gene_in,rv$obj)
    if(!length(gns)){showNotification("No genes found.",type="warning");return()}
    obj<-apply_ident(rv$obj,input$gene_ident%||%"cluster")
    if(!is.null(input$gene_clust_filter) && length(input$gene_clust_filter)) obj <- subset(obj, idents=input$gene_clust_filter)
    if(input$gene_type!="dot") obj <- filter_by_samples(obj, input$active_samples)
    id_labels<-unique(as.character(Idents(obj)))
    if(input$gene_type=="violin"&&length(gns)>1){
      n_id <- length(id_labels)
      w_scaled <- max(5.5, n_id*0.9, max_label_len(id_labels)*0.15)
      ht<-max(500,length(gns)*200)+max(90,max_label_len(id_labels)*6)
      dims<-list(w=w_scaled,h=max(6,length(gns)*2.8)+(max(90,max_label_len(id_labels)*6)/72))
    } else if(input$gene_type=="dot"){
      n_id<-length(id_labels)
      # Same sparse-case floor shrink as dotplot_dims (see its comment,
      # including why sqrt rather than a linear ramp), calibrated separately
      # for this tab's own 400px/10in/6in floors and its own dot.scale
      # (6, fixed -- untouched here).
      ht<-max(400*sqrt(pmin(1,n_id/4)),n_id*28,max_label_len(id_labels)*4)
      dims<-list(w=max(10*sqrt(pmin(1,length(gns)/4)),length(gns)*0.9),
                 h=max(6*sqrt(pmin(1,n_id/4)),n_id*0.4))
    } else if(input$gene_type=="feature"){
      ncf<-min(3,length(gns)); nrf<-ceiling(length(gns)/ncf)
      ht<-max(400,nrf*320); dims<-list(w=max(8,ncf*4),h=max(5,nrf*4))
    } else {
      ht<-560; dims<-list(w=12,h=8)
    }
    rv$cur_gene_dims<-dims
    output$gene_p_ui<-renderUI(
      div(style="overflow:auto; max-width:100%; max-height:80vh;",
          plotOutput("gene_p_inner",height=paste0(ht,"px"),width=paste0(px_from_in(dims$w),"px")))
    )
    output$gene_p_inner<-renderPlot({
      p<-if(input$gene_type=="feature"){
        lims<-umap_full_xylim(rv$obj)
        FeaturePlot(obj,features=gns,reduction="umap",pt.size=0.4,ncol=min(3,length(gns))) &  # explicit reduction
          coord_cartesian(xlim=lims$x, ylim=lims$y)
      }
      else if(input$gene_type=="violin")
        if(length(gns)==1) VlnPlot(obj,features=gns,pt.size=0.05,fill.by="ident") &
          theme(plot.margin=margin(t=60,r=20,b=max(90,max_label_len(id_labels)*6),l=40,unit="pt")) &
          fixed_ident_fill_scale(Idents(obj))
        else VlnPlot(obj,features=gns,pt.size=0.05,stack=TRUE,flip=TRUE,fill.by="ident") &
          theme(plot.margin=margin(t=60,r=max(20, max_label_len(gns)*4),b=max(90,max_label_len(id_labels)*6),l=40,unit="pt"), panel.spacing=unit(0.6,"cm")) &
          fixed_ident_fill_scale(Idents(obj))
      else DotPlot(obj,features=gns)+theme(axis.text.x=element_text(angle=45,hjust=1),
                                            plot.margin=margin(t=10,r=20,b=max(60,max_label_len(id_labels)*3),l=20,unit="pt"))
      p<-p&min_text_theme
      if(input$gene_type=="violin") p<-p&scale_y_continuous(expand=expansion(mult=c(0.05,0.15)))
      if(input$gene_type!="feature") p<-p&theme(axis.text.x=element_text(size=if(length(id_labels)<=28) 13 else pmax(9, 13-(length(id_labels)-28)*0.3)), axis.text.y=element_text(size=if(length(id_labels)<=28) 13 else pmax(9, 13-(length(id_labels)-28)*0.3)), axis.title=element_text(size=12))
      rv$cur_genes<-gns; rv$cur_gene_plot<-p; p
    })
  })

  # =========================================================================
  # OUTPUT FLAGS
  # =========================================================================
  output$data_loaded  <- reactive(isTRUE(rv$data_loaded))
  output$obj_raw_avail <- reactive(!is.null(rv$obj_raw))
  output$is_multi     <- reactive(isTRUE(rv$is_multi))
  output$qc_applied   <- reactive(isTRUE(rv$qc_applied))
  output$integ_done   <- reactive(isTRUE(rv$integ_done))
  output$has_integ_preview <- reactive(!is.null(rv$obj) && "umap.integrated.preview"%in%names(rv$obj@reductions))
  output$clust_done   <- reactive(isTRUE(rv$clust_done))
  output$markers_done <- reactive(isTRUE(rv$markers_done))
  output$sctype_done  <- reactive(isTRUE(rv$sctype_done))
  output$manual_annot_done <- reactive(isTRUE(rv$manual_annot_done))
  output$sig_done     <- reactive(isTRUE(rv$sig_done))
  output$cons_done    <- reactive(isTRUE(rv$cons_done))
  output$pb_done      <- reactive(isTRUE(rv$pb_done))
  output$de_done      <- reactive(isTRUE(rv$de_done))
  output$gsea_done    <- reactive(isTRUE(rv$gsea_done))
  output$ora_done     <- reactive(isTRUE(rv$ora_done))
  output$subprep_done  <- reactive(isTRUE(rv$subprep_done))
  output$sub_integ_done <- reactive(isTRUE(rv$sub_integ_done))
  output$subclust_done <- reactive(isTRUE(rv$subclust_done))
  output$sub_multi_parent <- reactive(isTRUE(length(rv$sub_parent_clusters)>1))
  output$sub_markers_done <- reactive(isTRUE(rv$sub_markers_done))
  for(x in c("data_loaded","obj_raw_avail","is_multi","qc_applied","integ_done","has_integ_preview","clust_done",
              "markers_done","sctype_done","manual_annot_done","sig_done","cons_done","pb_done","de_done","gsea_done","ora_done",
              "subprep_done","sub_integ_done","subclust_done","sub_multi_parent","sub_markers_done"))
    outputOptions(output,x,suspendWhenHidden=FALSE)

  # =========================================================================
  # DYNAMIC UIs
  # =========================================================================
  get_ident_ch <- function(){
    ch<-c("Cluster Number"="cluster")
    if(isTRUE(rv$sctype_done))       ch<-c(ch,"ScType Label"="sctype")
    if(isTRUE(rv$manual_annot_done)) ch<-c(ch,"Manual Label"="manual")
    if(!is.null(rv$obj)&&"subcluster_labels"%in%colnames(rv$obj@meta.data)) ch<-c(ch,"Subcluster Label"="subcluster")
    ch
  }
  output$clust_ident_ui <- renderUI({req(rv$clust_done);selectInput("clust_ident","Color identity:",get_ident_ch())})
  output$gene_ident_ui  <- renderUI({req(rv$clust_done);selectInput("gene_ident","Group by:",get_ident_ch())})
  output$gene_clust_filter_ui <- renderUI({
    req(rv$clust_done)
    ch <- nat_sort(unique(as.character(rv$obj@meta.data[[get_id_col(input$gene_ident%||%"cluster")]])))
    selectizeInput("gene_clust_filter","Clusters to include (optional):",choices=ch,multiple=TRUE)
  })
  output$mkr_group_by_txt   <- renderText(paste("Grouping by:", names(get_ident_ch())[get_ident_ch()==(input$clust_ident%||%"cluster")]))
  output$sig_group_by_txt   <- renderText(paste("Grouping by:", names(get_ident_ch())[get_ident_ch()==(input$clust_ident%||%"cluster")]))
  output$mannot_group_by_txt <- renderText(paste("Grouping by:", names(get_ident_ch())[get_ident_ch()==(input$clust_ident%||%"cluster")]))
  observe({req(rv$clust_done);updateSelectInput(session,"comp_ident",choices=get_ident_ch())})
  observe({req(rv$is_multi,length(rv$samp_names));updateSelectInput(session,"pre_split_sample",choices=rv$samp_names,selected=rv$samp_names[1])})
  observe({req(rv$is_multi,length(rv$samp_names));updateSelectInput(session,"cu_split_sample",choices=rv$samp_names,selected=rv$samp_names[1])})

  get_s1 <- function() if(!is.null(input$comp_s1)&&nchar(input$comp_s1)) input$comp_s1 else rv$samp_names[1]
  get_s2 <- function() if(!is.null(input$comp_s2)&&nchar(input$comp_s2)) input$comp_s2 else rv$samp_names[min(2,length(rv$samp_names))]

  output$de_samp_ui <- renderUI({
    req(rv$is_multi); sn<-rv$samp_names
    if(length(sn)<=2)
      div(class="ib",HTML(paste0("Comparing: <b>",sn[1],"</b> (Numerator / ident.1) vs <b>",sn[2],"</b> (Denominator / ident.2). avg_log2FC = log2(",sn[1]," / ",sn[2],"). Positive = higher in ",sn[1],".")))
    else tagList(
      div(class="ib","Select comparison pair. avg_log2FC = log2(Numerator / Denominator)."),
      fluidRow(
        column(4,selectInput("comp_s1","Numerator — ident.1 (A):",sn,selected=sn[1])),
        column(4,selectInput("comp_s2","Denominator — ident.2 (B):",sn,selected=sn[2])),
        column(4,uiOutput("de_formula_live"))
      )
    )
  })
  output$de_formula_live <- renderUI({div(class="ob",paste0("log2(",get_s1()," / ",get_s2(),")"))})
  output$de_formula_ui   <- renderUI({s1<-get_s1();s2<-get_s2();div(class="ib",HTML(paste0("avg_log2FC = log2(<b>",s1,"</b> / <b>",s2,"</b>). Positive = higher in ",s1,". Negative = higher in ",s2,".")))})
  output$mat_samp_ui     <- renderUI({req(rv$is_multi);selectInput("mat_samp","Sample filter:",c("All"="all",rv$samp_names),selected="all")})
  output$clust_reduc_label <- renderUI({req(rv$integ_done);paste0("Integrated UMAP | reduction: ",rv$reduc_use," | method: ",rv$integ_method)})
  output$comp_reduc_label  <- renderUI({req(rv$clust_done);paste0("Using integrated UMAP (reduction: ",rv$reduc_use,"). Expression = log-normalized counts.")})
  output$cons_clust_ui <- renderUI({req(rv$clust_done);selectInput("cons_clust","Cluster:",nat_sort(unique(as.character(rv$obj@meta.data[[get_id_col(input$comp_ident)]]))))})
  output$pb_ct_ui      <- renderUI({req(rv$clust_done);ch<-nat_sort(unique(as.character(rv$obj@meta.data[[get_id_col(input$comp_ident)]])));selectInput("pb_ct","Cell type / cluster:",ch,selected=ch[1],multiple=FALSE)})
  output$de_clust_ui   <- renderUI({req(rv$clust_done);ch<-nat_sort(unique(as.character(rv$obj@meta.data[[get_id_col(input$comp_ident)]])));selectInput("de_clust","Cluster / cell type:",ch)})
  # dev/diagnostic only — confirms cache entries accumulate instead of being overwritten; remove when no longer needed
  output$de_cache_debug <- renderPrint({ names(rv$de_cache) })
  output$mkr_clust_sel_ui <- renderUI({req(rv$markers_done);selectInput("mkr_clust_sel","Cluster:",nat_sort(unique(as.character(rv$markers$cluster))))})

  # =========================================================================
  # RENDER OUTPUTS
  # =========================================================================
  get_disp_obj <- function() rv$obj_joined %||% rv$obj

  output$data_summary <- renderPrint({
    req(get_disp_obj()); obj<-get_disp_obj()
    cat("=== Dataset Overview ===\n")
    cat("Cells:",ncol(obj),"| Genes:",nrow(obj),"\n")
    if(rv$is_multi){
      cat("\nSamples:\n"); t<-table(obj$sample)
      if(isTRUE(rv$qc_applied)&&!is.null(rv$obj_raw)){
        t_raw<-table(rv$obj_raw$sample)
        for(n in names(t)){
          kept<-t[[n]]
          if(n%in%names(t_raw)){
            orig<-t_raw[[n]]; removed<-orig-kept; pct<-round(100*kept/orig,1)
            cat(" ",n,": ",kept," cells kept (removed ",removed,", ",pct,"% retained)\n",sep="")
          } else cat(" ",n,": ",kept," cells\n",sep="")
        }
      } else if(is.null(rv$obj_raw)){
        for(n in names(t)) cat(" ",n,": ",t[[n]]," cells (loaded from saved object, pre-filter counts not available)\n",sep="")
      } else {
        for(n in names(t)) cat(" ",n,":",t[[n]],"cells\n")
      }
    }
    cat("\nMedian features/cell:",median(obj$nFeature_RNA),"\n")
    cat("Median UMIs/cell:",median(obj$nCount_RNA),"\n")
    cat("Median %MT:",round(median(obj$percent.mt),2),"%\n")
  })
  output$meta_dt <- renderDT({req(get_disp_obj());datatable(get_disp_obj()@meta.data,options=list(pageLength=25,scrollX=TRUE),rownames=TRUE)})
  output$samp_dist <- renderPlot({req(rv$obj,rv$is_multi)
    df<-as.data.frame(table(Sample=get_disp_obj()$sample))
    ggplot(df,aes(Sample,Freq,fill=Sample))+geom_col(width=0.6)+geom_text(aes(label=Freq),vjust=-0.3,size=4)+
      theme_minimal(base_size=13)+labs(y="Cells",x=NULL)+theme(legend.position="none",axis.text.x=element_text(angle=30,hjust=1))})

  # Counts matrix — LayerData for safe access
  get_mat <- function(obj,type) {
    tryCatch(switch(type,
      "counts"=as.matrix(LayerData(obj,layer="counts",assay="RNA")),
      "data"  =as.matrix(LayerData(obj,layer="data",assay="RNA")),
      "scale.data"=as.matrix(LayerData(obj,layer="scale.data",assay="RNA"))),
    error=function(e){showNotification(paste("Layer not available:",e$message),type="warning");NULL})
  }
  get_mat_obj <- function(){
    obj<-rv$obj
    if(rv$is_multi&&!is.null(input$mat_samp)&&input$mat_samp!="all")
      obj<-subset(obj,cells=colnames(obj)[obj$sample==input$mat_samp])
    obj
  }
  render_mat_dt <- function(){
    req(rv$obj); mat<-get_mat(get_mat_obj(),input$mat_type%||%"counts"); req(!is.null(mat))
    ng<-min(input$mat_genes%||%100,nrow(mat)); nc<-min(input$mat_cells%||%100,ncol(mat))
    df<-cbind(Gene=rownames(mat)[1:ng],as.data.frame(round(mat[1:ng,1:nc],3)))
    datatable(df,options=list(pageLength=20,scrollX=TRUE,scrollY="480px"),rownames=FALSE)
  }
  output$mat_dt <- renderDT({req(rv$obj);render_mat_dt()})
  observeEvent(input$update_mat,{output$mat_dt<-renderDT(render_mat_dt())})

  get_qc_obj <- function(){
    if(isTRUE(rv$qc_applied)&&!is.null(input$qc_view)&&input$qc_view=="after"&&!is.null(rv$obj_filt)) rv$obj_filt else rv$obj_raw
  }
  qc_missing_msg <- "QC plots not available when loading from a saved .rds file — the object was already filtered before saving."
  qc_unavailable_plot <- function() ggplot()+theme_void()+xlim(-1,1)+ylim(-1,1)+annotate("text",x=0,y=0,label=qc_missing_msg,size=5)
  output$qc_vln  <- renderPlot({req(rv$data_loaded);if(is.null(rv$obj_raw)) return(qc_unavailable_plot());VlnPlot(get_qc_obj(),features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,pt.size=0.05)})
  output$qc_scat <- renderPlot({req(rv$data_loaded);if(is.null(rv$obj_raw)) return(qc_unavailable_plot());FeatureScatter(get_qc_obj(),"nCount_RNA","percent.mt")+FeatureScatter(get_qc_obj(),"nCount_RNA","nFeature_RNA")})
  output$qc_samp <- renderPlot({req(rv$data_loaded,rv$is_multi);if(is.null(rv$obj_raw)) return(qc_unavailable_plot());VlnPlot(get_qc_obj(),features=c("nFeature_RNA","nCount_RNA","percent.mt"),group.by="sample",ncol=3,pt.size=0.05)})
  output$filt_summ <- renderPrint({req(rv$data_loaded);if(is.null(rv$obj_raw)){cat(qc_missing_msg,"\n");return(invisible(NULL))};req(rv$obj_filt);b<-ncol(rv$obj_raw);a<-ncol(rv$obj_filt);cat("Before:",b,"| After:",a,"| Removed:",b-a,"(",round((b-a)/b*100,1),"%)\n")})
  output$pca_p    <- renderPlot({req(rv$obj,rv$qc_applied);DimPlot(rv$obj,reduction="pca")+NoLegend()+ggtitle("PCA")})
  output$pca_samp <- renderPlot({req(rv$obj,rv$qc_applied,rv$is_multi);obj<-filter_by_samples(rv$obj,input$active_samples);DimPlot(obj,reduction="pca",group.by="sample")+ggtitle("PCA by Sample")})
  output$elbow    <- renderPlot({req(rv$obj,rv$qc_applied);ElbowPlot(rv$obj,ndims=40)})
  output$pca_load <- renderPlot({req(rv$obj,rv$qc_applied);VizDimLoadings(rv$obj,dims=1:min(input$n_load%||%4,ncol(rv$obj[["pca"]])),reduction="pca",ncol=2)})
  output$pca_heat <- renderPlot({req(rv$obj,rv$qc_applied);DimHeatmap(rv$obj,dims=1:min(input$n_heat%||%12,ncol(rv$obj[["pca"]])),cells=500,balanced=TRUE,ncol=3)})
  output$integ_summ <- renderPrint({req(rv$integ_done);cat("Method:",rv$integ_method,"| Reduction:",rv$reduc_use,"\n");cat("Available:",paste(names(rv$obj@reductions),collapse=", "),"\n")})
  integ_view_reduction <- function(){
    if(!is.null(input$integ_view) && input$integ_view=="post" && "umap.integrated.preview"%in%names(rv$obj@reductions)) "umap.integrated.preview" else "umap.unintegrated"
  }
  integ_view_title <- function() if(identical(integ_view_reduction(),"umap.integrated.preview")) "Post-integration (preview)" else "Pre-integration"
  output$pre_color_ui <- renderUI({
    req(rv$obj,rv$integ_done,rv$is_multi)
    dd <- repel_dims(unique(as.character(rv$obj$sample)))
    div(style="overflow:auto; max-width:100%; max-height:80vh;",
        plotOutput("pre_color", height=paste0(px_from_in(dd$h),"px"), width=paste0(px_from_in(dd$w),"px")))
  })
  output$pre_split_ui <- renderUI({
    req(rv$obj,rv$integ_done,rv$is_multi)
    dd <- repel_dims(unique(as.character(rv$obj$sample)))
    div(style="overflow:auto; max-width:100%; max-height:80vh;",
        plotOutput("pre_split", height=paste0(px_from_in(dd$h),"px"), width=paste0(px_from_in(dd$w),"px")))
  })
  output$pre_color  <- renderPlot({
    req(rv$obj,rv$integ_done,rv$is_multi,"umap.unintegrated"%in%names(rv$obj@reductions))
    obj <- filter_by_samples(rv$obj, input$active_samples)
    p <- umap_plot(obj,group.by="sample",title=integ_view_title(),pt.size=0.3,reduction=integ_view_reduction()) +
         fixed_ident_scale(rv$obj$sample)
    grid::grid.draw(fix_panel_size(p))
  })
  output$pre_split  <- renderPlot({
    req(rv$obj,rv$integ_done,rv$is_multi,input$pre_split_sample,"umap.unintegrated"%in%names(rv$obj@reductions))
    obj_s <- subset(rv$obj, subset = sample==input$pre_split_sample)
    obj_s <- filter_by_samples(obj_s, input$active_samples)
    p <- umap_plot(obj_s,group.by="sample",title=paste0(integ_view_title()," — ",input$pre_split_sample),pt.size=0.3,reduction=integ_view_reduction()) +
         fixed_ident_scale(rv$obj$sample)
    grid::grid.draw(fix_panel_size(p))
  })
  output$umap_p_ui <- renderUI({req(rv$obj,rv$clust_done)
    obj<-apply_ident(rv$obj,input$clust_ident%||%"cluster")
    dd<-repel_dims(unique(as.character(Idents(obj))))
    plotOutput("umap_p",height=paste0(px_from_in(dd$h),"px"),width=paste0(px_from_in(dd$w),"px"))})
  output$umap_p <- renderPlot({
    req(rv$obj,rv$clust_done)
    obj<-apply_ident(rv$obj,input$clust_ident%||%"cluster")
    obj<-filter_by_samples(obj,input$active_samples)
    if(!is.null(input$umap_col)&&input$umap_col=="sample")
      p <- umap_plot(obj,group.by="sample",title="UMAP by Sample",pt.size=0.4) + fixed_ident_scale(rv$obj$sample)
    else
      p <- umap_plot(obj,label=TRUE,pt.size=0.4,title=paste("UMAP | Res:",input$resolution,"| Method:",rv$integ_method)) +
           fixed_ident_scale(rv$obj@meta.data[[get_id_col(input$clust_ident%||%"cluster")]])
    grid::grid.draw(fix_panel_size(p))
  })
  output$clust_stats_dt <- renderDT({req(rv$obj,rv$clust_done)
    id_col <- get_id_col(input$clust_ident%||%"cluster")
    t<-table(rv$obj@meta.data[[id_col]])
    t<-t[nat_sort(names(t))]
    df<-data.frame(Cluster=names(t),Count=as.numeric(t),Pct=round(as.numeric(t)/sum(t)*100,2))
    if(rv$is_multi) for(s in rv$samp_names){cs<-rv$obj@meta.data[rv$obj$sample==s,id_col];ts<-table(cs);df[[paste0("n_",s)]]<-as.numeric(ts[df$Cluster])}
    datatable(df,options=list(pageLength=30))})
  output$cqc_vln  <- renderPlot({req(rv$obj,rv$clust_done);VlnPlot(rv$obj,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,pt.size=0.05)})
  output$cqc_scat <- renderPlot({req(rv$obj,rv$clust_done);FeatureScatter(rv$obj,"nCount_RNA","percent.mt")+FeatureScatter(rv$obj,"nCount_RNA","nFeature_RNA")})
  output$mkr_heat_ui <- renderUI({
    req(rv$markers_done, rv$markers)
    top <- rv$markers %>% group_by(cluster) %>% top_n(n=input$top_n %||% 10, wt=avg_log2FC)
    ng  <- length(unique(top$gene))
    plotOutput("mkr_heat", height=paste0(max(500, ng*18),"px"))
  })
  output$mkr_heat <- renderPlot({req(rv$markers,rv$obj);top<-rv$markers%>%group_by(cluster)%>%top_n(n=input$top_n,wt=avg_log2FC);nc<-length(unique(top$cluster));ng<-length(unique(top$gene));build_mkr_heatmap(rv$obj,top$gene,nc,ng,group_col=get_id_col(input$clust_ident%||%"cluster"))})
  output$mkr_top_dt <- renderDT({req(rv$markers);top<-rv$markers%>%group_by(cluster)%>%top_n(n=input$top_n,wt=avg_log2FC)%>%arrange(cluster,desc(avg_log2FC));datatable(top,filter="top",options=list(pageLength=25,scrollX=TRUE))})
  output$mkr_cluster_dt <- renderDT({req(rv$markers,input$mkr_clust_sel)
    d<-rv$markers%>%filter(cluster==input$mkr_clust_sel)%>%arrange(desc(avg_log2FC))
    if(!is.null(input$mkr_show_n)&&input$mkr_show_n>0) d<-head(d,input$mkr_show_n)
    datatable(d,filter="top",options=list(pageLength=25,scrollX=TRUE))})
  output$sub_mkr_heat_ui <- renderUI({
    req(rv$sub_markers_done, rv$sub_markers)
    top <- rv$sub_markers %>% group_by(cluster) %>% top_n(n=input$sub_top_n %||% 10, wt=avg_log2FC)
    ng  <- length(unique(top$gene))
    nc  <- length(unique(top$cluster))
    w_in <- max(10, nc*0.5 + max_label_len(as.character(unique(top$cluster)))*0.12) + 2
    plotOutput("sub_mkr_heat", height=paste0(max(500, ng*18),"px"), width=paste0(px_from_in(w_in),"px"))
  })
  output$sub_mkr_heat <- renderPlot({req(rv$sub_markers,rv$sub_obj);top<-rv$sub_markers%>%group_by(cluster)%>%top_n(n=input$sub_top_n,wt=avg_log2FC);nc<-length(unique(top$cluster));ng<-length(unique(top$gene));build_mkr_heatmap(rv$sub_obj,top$gene,nc,ng,group_col="subcluster_preview_label")}, res=72)
  output$sub_mkr_top_dt <- renderDT({req(rv$sub_markers);top<-rv$sub_markers%>%group_by(cluster)%>%top_n(n=input$sub_top_n,wt=avg_log2FC)%>%arrange(cluster,desc(avg_log2FC));datatable(top,filter="top",options=list(pageLength=25,scrollX=TRUE))})
  output$sctype_umap_ui <- renderUI({req(rv$obj,rv$sctype_done)
    dd<-repel_dims(unique(as.character(rv$obj$sctype_classification)))
    plotOutput("sctype_umap",height=paste0(px_from_in(dd$h),"px"),width=paste0(px_from_in(dd$w),"px"))})
  output$sctype_umap <- renderPlot({
    req(rv$obj,rv$sctype_done)
    obj <- filter_by_samples(rv$obj, input$active_samples)
    p <- umap_plot(obj,group.by="sctype_classification",label=TRUE,pt.size=0.4,title="ScType Annotation") +
         fixed_ident_scale(rv$obj$sctype_classification)
    grid::grid.draw(fix_panel_size(p))
  })
  output$sctype_tbl <- renderDT({ req(rv$sctype_res)
    datatable(rv$sctype_res %>%
                select(cluster, type, scores, ncells, marker_genes) %>%
                mutate(scores=round(scores, 2)),
              options=list(pageLength=30),
              rownames=FALSE) })
  output$manual_annot_umap_ui <- renderUI({req(rv$obj,rv$manual_annot_done)
    dd<-repel_dims(unique(as.character(rv$obj$manual_annotation)))
    plotOutput("manual_annot_umap",height=paste0(px_from_in(dd$h),"px"),width=paste0(px_from_in(dd$w),"px"))})
  output$manual_annot_umap <- renderPlot({
    req(rv$obj,rv$manual_annot_done)
    obj <- filter_by_samples(rv$obj, input$active_samples)
    p <- umap_plot(obj,group.by="manual_annotation",label=TRUE,pt.size=0.4,title="Manual Annotation") +
         fixed_ident_scale(rv$obj$manual_annotation)
    grid::grid.draw(fix_panel_size(p))
  })
  output$cu_samp_ui <- renderUI({req(rv$obj,rv$clust_done,rv$is_multi)
    dd<-repel_dims(rv$samp_names)
    plotOutput("cu_samp",height=paste0(px_from_in(dd$h),"px"),width=paste0(px_from_in(dd$w),"px"))})
  output$cu_samp  <- renderPlot({
    req(rv$obj,rv$clust_done,rv$is_multi)
    obj <- filter_by_samples(rv$obj, input$active_samples)
    p <- umap_plot(obj,group.by="sample",pt.size=0.3,title="UMAP by Sample") + fixed_ident_scale(rv$obj$sample)
    grid::grid.draw(fix_panel_size(p))
  })
  output$cu_split_ui <- renderUI({req(rv$obj,rv$clust_done,rv$is_multi)
    id_labels<-unique(as.character(rv$obj@meta.data[[get_id_col(input$comp_ident)]]))
    dd<-repel_dims(id_labels)
    plotOutput("cu_split",height=paste0(px_from_in(dd$h),"px"),width=paste0(px_from_in(dd$w),"px"))})
  output$cu_split <- renderPlot({
    req(rv$obj,rv$clust_done,rv$is_multi,input$cu_split_sample)
    color_by_sample <- identical(input$cu_split_color_by%||%"identity","sample")
    obj_s <- subset(apply_ident(rv$obj,input$comp_ident%||%"cluster"), subset = sample==input$cu_split_sample)
    obj_s <- filter_by_samples(obj_s, input$active_samples)
    if(color_by_sample){
      p <- umap_plot(obj_s, group.by="sample", pt.size=0.3, title=paste0("UMAP by Sample — ",input$cu_split_sample)) +
           fixed_ident_scale(rv$obj$sample)
    } else {
      p <- umap_plot(obj_s, label=TRUE, pt.size=0.3, title=paste0("UMAP by Identity — ",input$cu_split_sample)) +
           fixed_ident_scale(rv$obj@meta.data[[get_id_col(input$comp_ident%||%"cluster")]])
    }
    grid::grid.draw(fix_panel_size(p))
  })
  output$cu_ident_ui <- renderUI({req(rv$obj,rv$clust_done,rv$is_multi)
    obj<-apply_ident(rv$obj,input$comp_ident%||%"cluster")
    dd<-repel_dims(unique(as.character(Idents(obj))))
    plotOutput("cu_ident",height=paste0(px_from_in(dd$h),"px"),width=paste0(px_from_in(dd$w),"px"))})
  output$cu_ident <- renderPlot({
    req(rv$obj,rv$clust_done,rv$is_multi)
    obj<-apply_ident(rv$obj,input$comp_ident%||%"cluster")
    obj<-filter_by_samples(obj,input$active_samples)
    p <- umap_plot(obj,label=TRUE,pt.size=0.3,title="UMAP by Identity") +
         fixed_ident_scale(rv$obj@meta.data[[get_id_col(input$comp_ident%||%"cluster")]])
    grid::grid.draw(fix_panel_size(p))
  })
  output$comp_bar_ui <- renderUI({req(rv$obj,rv$clust_done,rv$is_multi)
    id<-get_id_col(input$comp_ident); labs_v<-unique(as.character(rv$obj@meta.data[[id]]))
    dd<-barflip_dims(length(labs_v),labs_v)
    plotOutput("comp_bar",height=paste0(px_from_in(dd$h),"px"))})
  output$comp_bar <- renderPlot({req(rv$obj,rv$clust_done,rv$is_multi)
    id<-get_id_col(input$comp_ident)
    labs_v<-unique(as.character(rv$obj@meta.data[[id]]))
    df<-rv$obj@meta.data%>%group_by(identity=.data[[id]],sample)%>%summarise(n=n(),.groups="drop")%>%group_by(identity)%>%mutate(pct=n/sum(n)*100)
    ggplot(df,aes(reorder(identity,pct),pct,fill=sample))+geom_col(width=0.7)+coord_flip()+theme_minimal(base_size=13)+
      labs(title="Cluster Composition",x=NULL,y="% cells",fill="Sample")+
      theme(plot.margin=margin(t=10,r=10,b=10,l=max(10,max_label_len(labs_v)*1.5),unit="pt"))})
  output$comp_bar_persample_ui <- renderUI({req(rv$obj,rv$clust_done,rv$is_multi)
    id<-get_id_col(input$comp_ident); labs_v<-unique(as.character(rv$obj@meta.data[[id]]))
    dd<-bardodge_dims(length(labs_v),labs_v)
    plotOutput("comp_bar_persample",height=paste0(px_from_in(dd$h),"px"),width=paste0(px_from_in(dd$w),"px"))})
  output$comp_bar_persample <- renderPlot({
    req(rv$obj, rv$clust_done, rv$is_multi)
    id <- get_id_col(input$comp_ident)
    labs_v <- unique(as.character(rv$obj@meta.data[[id]]))
    df <- rv$obj@meta.data %>%
      group_by(sample, identity=.data[[id]]) %>%
      summarise(n=n(), .groups="drop") %>%
      group_by(sample) %>%
      mutate(pct = n/sum(n)*100)
    ggplot(df, aes(x=identity, y=pct, fill=sample)) +
      geom_col(position=position_dodge(width=0.8), width=0.75) +
      theme_minimal(base_size=13) +
      theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
            plot.margin=margin(t=10,r=10,b=max(40,max_label_len(labs_v)*2),l=10,unit="pt")) +
      labs(title="Cell Type Composition per Sample", x=NULL, y="% of sample cells", fill="Sample")
  })
  output$cons_dt <- renderDT({req(rv$cons_markers);datatable(rv$cons_markers,filter="top",options=list(pageLength=25,scrollX=TRUE))})
  output$pb_p <- renderPlot({
    req(rv$pb_obj, rv$pb_done, input$pb_ct)
    s1 <- get_s1(); s2 <- get_s2()
    ct <- as.character(input$pb_ct[1])
    all_ids <- colnames(rv$pb_obj)
    i1 <- find_pb_col(rv$pb_obj, ct, s1, rv$pb_id_col, rv$obj)
    i2 <- find_pb_col(rv$pb_obj, ct, s2, rv$pb_id_col, rv$obj)
    if(is.na(i1)||is.na(i2)){
      plot.new()
      text(0.5,0.5,paste0("Not found.\nAvailable:\n",paste(all_ids,collapse="\n")),cex=0.65)
      return()
    }
    hi <- if(nchar(trimws(input$pb_hi)) > 0)
            resolve_genes(input$pb_hi, rv$obj) else character(0)
    r_val <- tryCatch({
      pb_mat <- GetAssayData(rv$pb_obj, assay=DefaultAssay(rv$pb_obj), layer="data")
      if(i1 %in% colnames(pb_mat) && i2 %in% colnames(pb_mat))
        round(cor(as.numeric(pb_mat[,i1]), as.numeric(pb_mat[,i2]), use="complete.obs"), 3)
      else NA_real_
    }, error=function(e) NA_real_)
    p <- CellScatter(rv$pb_obj, i1, i2,
                     highlight = if(length(hi)) hi else NULL)
    if(length(hi)) p <- LabelPoints(plot=p, points=hi, repel=TRUE)
    p + labs(title=paste0("Pseudobulk: ",ct,"  |  ",s1," vs ",s2),
             subtitle=if(!is.na(r_val)) paste0("Pearson r = ",r_val) else NULL) +
      theme(plot.subtitle=element_text(size=13, face="bold", color="#2c3e50", hjust=0.5))
  })
  output$pb_warn_ui <- renderUI({
    req(rv$pb_done)
    if(!is.null(rv$pb_id_col) && rv$pb_id_col != get_id_col(input$comp_ident))
      div(class="wb", paste0("Warning: pseudobulk was built with grouping '",rv$pb_id_col,
                             "' but current identity is '",get_id_col(input$comp_ident),
                             "'. Re-run Pseudobulk for accurate results."))
    else NULL
  })
  output$de_dt <- renderDT({
    req(rv$de_res)
    lfc_t  <- as.numeric(input$volc_lfc)  %||% 1
    pval_t <- as.numeric(input$volc_pval) %||% 0.05
    df <- as.data.frame(rv$de_res)
    df$Significant <- ifelse(abs(df$avg_log2FC) >= lfc_t & df$p_val_adj <= pval_t, "Yes", "No")
    datatable(df, filter="top", options=list(pageLength=25, scrollX=TRUE)) %>%
      formatStyle("Significant",
        backgroundColor=styleEqual(c("Yes","No"), c("#d4edda","#f5f5f5")),
        color=styleEqual(c("Yes","No"), c("#155724","#888888")),
        fontWeight=styleEqual(c("Yes","No"), c("bold","normal")))
  })
  output$volc_p <- renderPlot({
    req(rv$de_res, rv$de_done, nrow(rv$de_res) > 0)
    tryCatch(make_volcano(rv$de_res,lfc_thresh=as.numeric(input$volc_lfc)%||%1,pval_thresh=as.numeric(input$volc_pval)%||%0.05,
                          top_n=as.integer(input$volc_n)%||%10,
                          title=paste0("DE: ",input$de_clust%||%""," | log\u2082(",rv$de_s1," / ",rv$de_s2,")"),
                          s1=rv$de_s1,s2=rv$de_s2),
             error=function(e){plot.new();text(0.5,0.5,paste("Volcano error:\n",e$message),cex=0.9,col="red")})
  })

  # ORA renders — shared helper combines both directions into one data frame
  # (direction column added), mirroring gsea_top_df's pattern so the in-app
  # render and both downloads always stay in sync.
  ora_combined <- function(){
    req(rv$ora_done)
    lab_up <- paste0("Higher in ", rv$ora_s1); lab_dn <- paste0("Higher in ", rv$ora_s2)
    parts <- list()
    if(!is.null(rv$ora_res_up)&&nrow(rv$ora_res_up)>0)
      parts[[length(parts)+1]] <- cbind(rv$ora_res_up, direction=lab_up, stringsAsFactors=FALSE)
    if(!is.null(rv$ora_res_dn)&&nrow(rv$ora_res_dn)>0)
      parts[[length(parts)+1]] <- cbind(rv$ora_res_dn, direction=lab_dn, stringsAsFactors=FALSE)
    if(!length(parts)) return(NULL)
    do.call(rbind, parts)
  }
  output$ora_has_results <- reactive(!is.null(ora_combined()))
  outputOptions(output,"ora_has_results",suspendWhenHidden=FALSE)
  output$ora_plot_ui <- renderUI({
    df <- ora_combined()
    if(is.null(df)) return(div(class="ib","No significant terms found."))
    dd <- barflip_dims(nrow(df), trunc_pw(df$Description))
    plotOutput("ora_plot", height=paste0(px_from_in(dd$h),"px"), width=paste0(px_from_in(dd$w),"px"))
  })
  output$ora_plot <- renderPlot({
    df <- ora_combined()
    req(!is.null(df))
    ora_make_plot(df, rv$ora_label, rv$ora_s1, rv$ora_s2, input$ora_plot_type %||% "bar")
  })

  # GSEA — shared helper so the in-app render and the download stay in sync
  gsea_top_df <- function(){
    req(rv$gsea_done,rv$gsea_res)
    topN<-max(1,as.integer(input$gsea_topn%||%15))
    fg_up<-rv$gsea_res%>%filter(padj<0.05,NES>0)%>%arrange(padj)%>%head(topN)
    fg_down<-rv$gsea_res%>%filter(padj<0.05,NES<0)%>%arrange(padj)%>%head(topN)
    bind_rows(fg_up,fg_down)
  }
  output$gsea_plot_ui <- renderUI({
    top_df<-gsea_top_df()
    if(!nrow(top_df)) return(plotOutput("gsea_plot",height="150px"))
    dd<-barflip_dims(nrow(top_df),trunc_pw(top_df$pathway))
    dd$w <- dd$w + gsea_legend_extra_in  # external patchwork legend column, see gsea_with_legend()
    plotOutput("gsea_plot",height=paste0(px_from_in(dd$h),"px"),width=paste0(px_from_in(dd$w),"px"))
  })
  output$gsea_plot <- renderPlot({
    top_df<-gsea_top_df()
    if(!nrow(top_df)) return(ggplot()+theme_void()+annotate("text",x=0,y=0,label="No significant pathways (padj<0.05) for either direction."))
    p<-gsea_make_plot(top_df,input$gsea_plot_type%||%"dot",rv$gsea_s1,rv$gsea_s2,rv$gsea_label)
    gsea_with_legend(p)
  })
  output$gsea_dt <- renderDT({
    req(rv$gsea_done,rv$gsea_res)
    df<-rv$gsea_res%>%mutate(leadingEdge=sapply(leadingEdge,function(x) paste(x,collapse=", ")),
                              pval=signif(pval,4),padj=signif(padj,4),ES=round(ES,4),NES=round(NES,4))
    datatable(df,filter="top",options=list(pageLength=25,scrollX=TRUE))
  })

  # ---- PATHWAY EXPLORER ----
  observeEvent(rv$gsea_res, {
    req(rv$gsea_done, rv$gsea_res)
    choices <- unique(rv$gsea_res$pathway)
    top_sig  <- head(choices[rv$gsea_res$padj < 0.05], 5)
    updateSelectizeInput(session, "pe_pathways", choices=choices,
                         selected=top_sig, server=TRUE)
  })

  pe_build_bar <- function(df, s1, s2, label) {
    n <- nrow(df)
    disp <- trunc_pw(df$pathway)
    df$pathway_disp <- disp
    df$sig_label <- ifelse(df$padj < 0.05, "*", "")
    fsz <- mkr_font(n, 9, 11, 200)
    nes_range <- max(abs(df$NES), na.rm=TRUE)
    ggplot(df, aes(x=reorder(pathway_disp, NES), y=NES, fill=NES > 0)) +
      geom_hline(yintercept=0, color="black", linewidth=0.8) +
      geom_col(alpha=0.85) +
      geom_text(aes(label=sig_label,
                    y=NES + sign(NES) * nes_range * 0.08),
                size=5, fontface="bold") +
      coord_flip() +
      scale_fill_manual(
        values=c("TRUE"="#e74c3c", "FALSE"="#3498db"),
        labels=c("TRUE"=paste("Higher in", s1),
                 "FALSE"=paste("Higher in", s2)),
        name="Direction") +
      labs(title=label, x=NULL, y="NES",
           subtitle=paste0("* padj<0.05  |  Positive NES = higher in ", s1)) +
      theme_minimal(base_size=12) +
      theme(axis.text.y=element_text(size=fsz,face="bold"),
            plot.margin=margin(t=10,r=10,b=10,l=max(80,max_label_len(disp)*5),
                               unit="pt"),
            plot.subtitle=element_text(size=10, color="#555")) +
      min_text_theme
  }

  output$pe_plot_ui <- renderUI({
    req(rv$gsea_done, rv$gsea_res, input$pe_pathways, length(input$pe_pathways) > 0)
    df <- rv$gsea_res %>% filter(pathway %in% input$pe_pathways)
    if (!nrow(df)) return(plotOutput("pe_plot", height="150px"))
    if (identical(input$pe_view, "trace")) {
      ht <- max(300, length(input$pe_pathways) * 220)
      return(plotOutput("pe_plot", height=paste0(ht, "px")))
    }
    dd <- barflip_dims(nrow(df), trunc_pw(df$pathway))
    plotOutput("pe_plot",
               height=paste0(px_from_in(dd$h), "px"),
               width =paste0(px_from_in(dd$w), "px"))
  })

  output$pe_plot <- renderPlot({
    req(rv$gsea_done, rv$gsea_res, input$pe_pathways, length(input$pe_pathways) > 0)
    df <- rv$gsea_res %>% filter(pathway %in% input$pe_pathways)
    if (!nrow(df)) {
      return(ggplot() + theme_void() +
               annotate("text", x=0, y=0, label="No matching pathways in current results."))
    }
    if (identical(input$pe_view %||% "barplot", "trace")) {
      req(!is.null(rv$gsea_pl), !is.null(rv$gsea_stats))
      plots <- lapply(input$pe_pathways, function(pw) {
        if (pw %in% names(rv$gsea_pl)) {
          fgsea::plotEnrichment(rv$gsea_pl[[pw]], rv$gsea_stats) +
            labs(title=trunc_pw(pw, max_len=100)) +
            theme_minimal(base_size=11) + min_text_theme
        } else NULL
      })
      plots <- Filter(Negate(is.null), plots)
      if (!length(plots))
        return(ggplot() + theme_void() +
                 annotate("text", x=0, y=0, label="Pathway not found in current gene set."))
      wrap_plots(plots, ncol=1)
    } else {
      pe_build_bar(df, rv$gsea_s1, rv$gsea_s2, rv$gsea_label)
    }
  })

  # =========================================================================
  # DOWNLOADS
  # =========================================================================
  output$dl_seurat <- downloadHandler(filename=function() paste0("scgui_obj_",Sys.Date(),".rds"),content=function(f) saveRDS(rv$obj,f))
  output$dl_meta_csv <- downloadHandler(filename=function() paste0("metadata_",Sys.Date(),".csv"),content=function(f) write.csv(as.data.frame(get_disp_obj()@meta.data),f))
  output$dl_meta_csv2 <- downloadHandler(filename=function() paste0("metadata_all_",Sys.Date(),".csv"),content=function(f) write.csv(as.data.frame(get_disp_obj()@meta.data),f))
  # TSV for matrix (avoids Excel column limit)
  output$dl_mat <- downloadHandler(
    filename=function() paste0("matrix_",input$mat_type%||%"counts","_",Sys.Date(),".tsv"),
    content=function(f){
      mat<-get_mat(get_mat_obj(),input$mat_type%||%"counts"); if(is.null(mat)) return()
      write.table(cbind(Gene=rownames(mat),as.data.frame(mat)),f,sep="\t",row.names=FALSE,quote=FALSE)})
  output$dl_mkr_top     <- downloadHandler(filename=function() paste0("top",input$top_n,"_markers_",Sys.Date(),".csv"),content=function(f) write.csv(rv$markers%>%group_by(cluster)%>%top_n(n=input$top_n,wt=avg_log2FC),f,row.names=FALSE))
  output$dl_mkr_all     <- downloadHandler(filename=function() paste0("all_markers_",Sys.Date(),".csv"),content=function(f) write.csv(rv$markers,f,row.names=FALSE))
  output$dl_mkr_cluster <- downloadHandler(filename=function() paste0("markers_cluster",input$mkr_clust_sel,"_",Sys.Date(),".csv"),content=function(f) write.csv(rv$markers%>%filter(cluster==input$mkr_clust_sel)%>%arrange(desc(avg_log2FC)),f,row.names=FALSE))
  output$dl_sub_mkr_top <- downloadHandler(filename=function() paste0("sub_top",input$sub_top_n,"_markers_",Sys.Date(),".csv"),content=function(f) write.csv(rv$sub_markers%>%group_by(cluster)%>%top_n(n=input$sub_top_n,wt=avg_log2FC),f,row.names=FALSE))
  output$dl_sub_mkr_all <- downloadHandler(filename=function() paste0("sub_all_markers_",Sys.Date(),".csv"),content=function(f) write.csv(rv$sub_markers,f,row.names=FALSE))
  output$dl_sctype_tbl  <- downloadHandler(filename=function() paste0("sctype_",Sys.Date(),".csv"),content=function(f) write.csv(as.data.frame(rv$sctype_res),f,row.names=FALSE))
  output$dl_comp_csv    <- downloadHandler(filename=function() paste0("composition_",Sys.Date(),".csv"),content=function(f){id<-get_id_col(input$comp_ident);df<-rv$obj@meta.data%>%group_by(identity=.data[[id]],sample)%>%summarise(n=n(),.groups="drop")%>%group_by(identity)%>%mutate(pct=round(n/sum(n)*100,2));write.csv(df,f,row.names=FALSE)})
  output$dl_cons        <- downloadHandler(filename=function() paste0("conserved_markers_",Sys.Date(),".csv"),content=function(f) write.csv(rv$cons_markers,f,row.names=FALSE))
  output$dl_de          <- downloadHandler(filename=function() paste0("DE_",input$de_clust%||%""    ,"_",rv$de_s1,"_vs_",rv$de_s2,"_",Sys.Date(),".csv"),content=function(f) write.csv(rv$de_res,f,row.names=FALSE))
  output$dl_clust_stats <- downloadHandler(filename=function() paste0("cluster_stats_",Sys.Date(),".csv"),content=function(f){id_col<-get_id_col(input$clust_ident%||%"cluster");t<-table(rv$obj@meta.data[[id_col]]);t<-t[nat_sort(names(t))];write.csv(data.frame(Cluster=names(t),Count=as.numeric(t),Pct=round(as.numeric(t)/sum(t)*100,2)),f,row.names=FALSE)})

  # 600 DPI plot downloads
  output$dl_samp_dist <- dl600(function(){df<-as.data.frame(table(Sample=get_disp_obj()$sample));ggplot(df,aes(Sample,Freq,fill=Sample))+geom_col(width=0.6)+geom_text(aes(label=Freq),vjust=-0.3)+theme_minimal()+labs(y="Cells",x=NULL)+theme(legend.position="none")},"sample_dist",8,5)
  output$dl_qc_vln    <- dl600(function() if(is.null(rv$obj_raw)) qc_unavailable_plot() else VlnPlot(get_qc_obj(),features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,pt.size=0.05),"qc_violin",14,6)
  output$dl_qc_scat   <- dl600(function() if(is.null(rv$obj_raw)) qc_unavailable_plot() else FeatureScatter(get_qc_obj(),"nCount_RNA","percent.mt")+FeatureScatter(get_qc_obj(),"nCount_RNA","nFeature_RNA"),"qc_scatter",12,5)
  output$dl_qc_samp   <- dl600(function() if(is.null(rv$obj_raw)) qc_unavailable_plot() else VlnPlot(get_qc_obj(),features=c("nFeature_RNA","nCount_RNA","percent.mt"),group.by="sample",ncol=3,pt.size=0.05),"qc_by_sample",14,6)
  output$dl_pca_p     <- dl600(function() DimPlot(rv$obj,reduction="pca")+NoLegend(),"pca")
  output$dl_pca_samp  <- dl600(function() {obj<-filter_by_samples(rv$obj,input$active_samples); DimPlot(obj,reduction="pca",group.by="sample")},"pca_sample")
  output$dl_elbow     <- dl600(function() ElbowPlot(rv$obj,ndims=40),"elbow",8,5)
  output$dl_pca_load  <- dl600(function() VizDimLoadings(rv$obj,dims=1:min(input$n_load%||%4,ncol(rv$obj[["pca"]])),reduction="pca",ncol=2),"pca_loadings",10,12)
  output$dl_pca_heat  <- downloadHandler(
    filename=function() paste0("pca_heatmap_",Sys.Date(),".png"),
    content=function(f){req(rv$obj,rv$qc_applied);n<-min(input$n_heat%||%12,ncol(rv$obj[["pca"]]));png(f,width=9000,height=12000,res=600);DimHeatmap(rv$obj,dims=1:n,cells=500,balanced=TRUE,ncol=3);dev.off()})
  output$dl_sub_elbow <- downloadHandler(
    filename=function() paste0("sub_elbow_",Sys.Date(),".png"),
    content=function(f){
      req(rv$sub_obj)
      ndims <- min(40,ncol(rv$sub_obj[["pca"]]))
      w_in <- max(8, ndims*0.15); h_in <- 5
      p <- ElbowPlot(rv$sub_obj,ndims=ndims)
      ggsave(f,p,width=w_in,height=h_in,dpi=600,bg="white")
    })
  output$dl_sub_pca_heat <- downloadHandler(
    filename=function() paste0("sub_pca_heatmap_",Sys.Date(),".png"),
    content=function(f){
      req(rv$sub_obj)
      n <- min(input$n_heat%||%12,ncol(rv$sub_obj[["pca"]]))
      nrow_grid <- ceiling(n/3)
      w_in <- 12; h_in <- max(6, nrow_grid*3)
      png(f,width=px_from_in(w_in,ppi=600),height=px_from_in(h_in,ppi=600),res=600)
      DimHeatmap(rv$sub_obj,dims=1:n,cells=500,balanced=TRUE,ncol=3)
      dev.off()
    })
  output$dl_sub_umap_pre <- dl600(function(){
    obj <- filter_by_samples(rv$sub_obj, input$active_samples)
    p <- umap_plot(obj,group.by="sample",pt.size=0.3,title=sub_integ_view_title(),reduction=sub_integ_view_reduction()) +
         fixed_ident_scale(rv$sub_obj$sample)
    fix_panel_size(p)
  },"sub_preint")
  output$dl_sub_umap <- downloadHandler(
    filename=function() paste0("subcluster_umap_",Sys.Date(),".png"),
    content=function(f){
      req(rv$sub_obj,rv$subclust_done)
      labs_v <- nat_sort(unique(as.character(Idents(rv$sub_obj))))
      dd <- repel_dims(labs_v)
      obj <- filter_by_samples(rv$sub_obj, input$active_samples)
      p <- umap_plot(obj,label=TRUE,pt.size=0.4,title="Subcluster UMAP") +
           fixed_ident_scale(Idents(rv$sub_obj))
      ggsave(f,fix_panel_size(p),width=dd$w,height=dd$h,dpi=600,bg="white",device="png")
    })
  output$dl_pre_color <- dl600(function() {obj<-filter_by_samples(rv$obj,input$active_samples); fix_panel_size(umap_plot(obj,group.by="sample",title=integ_view_title(),pt.size=0.3,reduction=integ_view_reduction()) + fixed_ident_scale(rv$obj$sample))},"preint_color")
  output$dl_pre_split <- dl600(function(){
    obj_s <- subset(rv$obj, subset = sample==input$pre_split_sample)
    obj_s <- filter_by_samples(obj_s, input$active_samples)
    p <- umap_plot(obj_s,group.by="sample",title=paste0(integ_view_title()," — ",input$pre_split_sample),pt.size=0.3,reduction=integ_view_reduction()) +
         fixed_ident_scale(rv$obj$sample)
    fix_panel_size(p)
  },"preint_split")
  output$dl_umap <- downloadHandler(
    filename=function() paste0("umap_",Sys.Date(),".png"),
    content=function(f){
      if(!is.null(input$umap_col) && input$umap_col=="sample"){
        dd <- repel_dims(unique(as.character(rv$obj$sample)))
        obj <- filter_by_samples(rv$obj, input$active_samples)
        p <- umap_plot(obj,group.by="sample",pt.size=0.4,title="UMAP by Sample")+
             guides(color=guide_legend(override.aes=list(size=4)))+
             fixed_ident_scale(rv$obj$sample)
        ggsave(f,fix_panel_size(p),width=dd$w,height=dd$h,dpi=600,bg="white",device="png")
      } else {
        obj <- apply_ident(rv$obj,input$clust_ident%||%"cluster")
        obj <- filter_by_samples(obj, input$active_samples)
        labs_v <- unique(as.character(Idents(obj)))
        dd  <- repel_dims(labs_v)
        p   <- umap_plot(obj,label=TRUE,pt.size=0.4)+
               guides(color=guide_legend(override.aes=list(size=4),ncol=if(length(labs_v)>28) 2 else 1))+
               fixed_ident_scale(rv$obj@meta.data[[get_id_col(input$clust_ident%||%"cluster")]])
        ggsave(f,fix_panel_size(p),width=dd$w,height=dd$h,dpi=600,bg="white",device="png")
      }
    })
  output$dl_cqc_vln   <- dl600(function() VlnPlot(rv$obj,features=c("nFeature_RNA","nCount_RNA","percent.mt"),ncol=3,pt.size=0.05),"cqc_violin",14,6)
  output$dl_cqc_scat  <- dl600(function() FeatureScatter(rv$obj,"nCount_RNA","percent.mt")+FeatureScatter(rv$obj,"nCount_RNA","nFeature_RNA"),"cqc_scatter")
  output$dl_gene <- downloadHandler(
    filename=function() paste0("gene_expr_",Sys.Date(),".png"),
    content=function(f){
      req(rv$cur_gene_plot)
      dims <- rv$cur_gene_dims %||% list(w=12,h=8)
      ggsave(f,rv$cur_gene_plot,width=min(dims$w,16),height=min(dims$h,14),
             dpi=600,bg="white",device="png")
    })
  output$dl_mkr_heat  <- downloadHandler(filename=function() paste0("marker_heatmap_",Sys.Date(),".png"),content=function(f){req(rv$markers,rv$obj);top<-rv$markers%>%group_by(cluster)%>%top_n(n=input$top_n,wt=avg_log2FC);nc<-length(unique(top$cluster));ng<-length(unique(top$gene));p<-build_mkr_heatmap(rv$obj,top$gene,nc,ng,group_col=get_id_col(input$clust_ident%||%"cluster"));ggsave(f,p,width=max(10,nc*0.9),height=min(20,max(8,ng*0.15)),dpi=600,limitsize=FALSE,bg="white")})
  output$dl_sub_mkr_heat <- downloadHandler(filename=function() paste0("sub_marker_heatmap_",Sys.Date(),".png"),content=function(f){req(rv$sub_markers,rv$sub_obj);top<-rv$sub_markers%>%group_by(cluster)%>%top_n(n=input$sub_top_n,wt=avg_log2FC);nc<-length(unique(top$cluster));ng<-length(unique(top$gene));p<-build_mkr_heatmap(rv$sub_obj,top$gene,nc,ng,group_col="subcluster_preview_label");ggsave(f,p,width=(max(10,nc*0.9+max_label_len(as.character(unique(top$cluster)))*0.12)+2)*(92/72),height=min(20,max(8,ng*0.15))*(92/72),dpi=600,limitsize=FALSE,bg="white")})
  output$dl_sctype_umap <- downloadHandler(filename=function() paste0("sctype_umap_",Sys.Date(),".png"),content=function(f){labs_v<-unique(as.character(rv$obj$sctype_classification));dd<-repel_dims(labs_v,base_w=12,base_h=9);obj<-filter_by_samples(rv$obj,input$active_samples);p<-umap_plot(obj,group.by="sctype_classification",label=TRUE,pt.size=0.4,title="ScType Annotation")+guides(color=guide_legend(override.aes=list(size=4),ncol=if(length(labs_v)>28) 2 else 1))+fixed_ident_scale(rv$obj$sctype_classification);ggsave(f,fix_panel_size(p),width=dd$w,height=dd$h,dpi=600,bg="white",device="png")})
  output$dl_manual_annot_umap <- downloadHandler(filename=function() paste0("manual_annot_umap_",Sys.Date(),".png"),content=function(f){labs_v<-unique(as.character(rv$obj$manual_annotation));dd<-repel_dims(labs_v,base_w=12,base_h=9);obj<-filter_by_samples(rv$obj,input$active_samples);p<-umap_plot(obj,group.by="manual_annotation",label=TRUE,pt.size=0.4,title="Manual Annotation")+guides(color=guide_legend(override.aes=list(size=4),ncol=if(length(labs_v)>28) 2 else 1))+fixed_ident_scale(rv$obj$manual_annotation);ggsave(f,fix_panel_size(p),width=dd$w,height=dd$h,dpi=600,bg="white",device="png")})
  output$dl_comp_bar_persample <- downloadHandler(
    filename=function() paste0("composition_persample_",Sys.Date(),".png"),
    content=function(f){
      req(rv$obj,rv$clust_done,rv$is_multi)
      id <- get_id_col(input$comp_ident)
      df <- rv$obj@meta.data %>%
        group_by(sample, identity=.data[[id]]) %>%
        summarise(n=n(), .groups="drop") %>%
        group_by(sample) %>%
        mutate(pct = n/sum(n)*100)
      n_clusters <- length(unique(df$identity))
      labs_v <- unique(as.character(df$identity))
      dd <- bardodge_dims(n_clusters, labs_v)
      p <- ggplot(df, aes(x=identity, y=pct, fill=sample)) +
        geom_col(position=position_dodge(width=0.8), width=0.75) +
        theme_minimal(base_size=13) +
        theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
              plot.margin=margin(t=10,r=10,b=max(40,max_label_len(labs_v)*2),l=10,unit="pt")) +
        labs(title="Cell Type Composition per Sample", x=NULL, y="% of sample cells", fill="Sample")
      ggsave(f, p, width=dd$w, height=dd$h, dpi=600, bg="white")
    })
  output$dl_comp_persample_csv <- downloadHandler(
    filename=function() paste0("composition_persample_",Sys.Date(),".csv"),
    content=function(f){
      req(rv$obj,rv$clust_done,rv$is_multi)
      id <- get_id_col(input$comp_ident)
      df <- rv$obj@meta.data %>%
        group_by(sample, identity=.data[[id]]) %>%
        summarise(n=n(), .groups="drop") %>%
        group_by(sample) %>%
        mutate(pct=round(n/sum(n)*100,2))
      write.csv(df, f, row.names=FALSE)
    })
  output$dl_sig_heat  <- downloadHandler(
    filename=function() paste0("signature_heatmap_",Sys.Date(),".png"),
    content=function(f){
      req(rv$sig_scores)
      mat <- rv$sig_scores
      zmat <- if(nrow(mat)>1) t(scale(t(mat))) else mat
      zmat[!is.finite(zmat)] <- 0
      col_margin <- max(8, max_label_len(colnames(zmat))*0.6)
      row_margin <- max(10, max_label_len(rownames(zmat))*0.6)
      w_in <- max(10, ncol(zmat)*0.45 + max_label_len(colnames(zmat))*0.12)
      h_in <- max(8, nrow(zmat)*0.45 + max_label_len(rownames(zmat))*0.12)
      png_scale <- 92/72  # matches on-screen's actual effective size: px_from_in() assumes
                           # 92 px/in for the on-screen container, but renderPlot renders at
                           # 72 dpi -- this ratio is the exact correction, not an eyeballed guess
      png(f, width=w_in*png_scale*600, height=h_in*png_scale*600, res=600)
      heatmap(zmat, scale="none",
              col=colorRampPalette(c("#3498db","white","#e74c3c"))(100),
              margins=c(col_margin, row_margin), main="Signature Scores (z-scored per signature)",
              cexRow=1.0, cexCol=1.0)
      dev.off()
    })
  output$dl_sig_csv   <- downloadHandler(
    filename=function() paste0("signature_scores_",Sys.Date(),".csv"),
    content=function(f){
      req(rv$sig_scores)
      df <- cbind(Signature=rownames(rv$sig_scores), as.data.frame(round(rv$sig_scores,4)))
      write.csv(df, f, row.names=FALSE)
    })
  output$dl_cu_samp <- downloadHandler(
    filename=function() paste0("comp_umap_sample_",Sys.Date(),".png"),
    content=function(f){
      dd <- repel_dims(rv$samp_names)
      obj <- filter_by_samples(rv$obj, input$active_samples)
      p <- umap_plot(obj,group.by="sample",pt.size=0.3,title="UMAP by Sample")+
           guides(color=guide_legend(override.aes=list(size=4)))+
           fixed_ident_scale(rv$obj$sample)
      ggsave(f,fix_panel_size(p),width=dd$w,height=dd$h,dpi=600,bg="white",device="png")
    })

  output$dl_cu_split <- downloadHandler(
    filename=function() paste0("comp_umap_split_",Sys.Date(),".png"),
    content=function(f){
      id_labels <- unique(as.character(rv$obj@meta.data[[get_id_col(input$comp_ident)]]))
      dd <- repel_dims(id_labels)
      color_by_sample <- identical(input$cu_split_color_by%||%"identity","sample")
      obj_s <- subset(apply_ident(rv$obj,input$comp_ident%||%"cluster"), subset = sample==input$cu_split_sample)
      obj_s <- filter_by_samples(obj_s, input$active_samples)
      if(color_by_sample){
        p <- umap_plot(obj_s, group.by="sample", pt.size=0.3, title=paste0("UMAP by Sample — ",input$cu_split_sample)) +
             fixed_ident_scale(rv$obj$sample)
      } else {
        p <- umap_plot(obj_s, label=TRUE, pt.size=0.3, title=paste0("UMAP by Identity — ",input$cu_split_sample)) +
             fixed_ident_scale(rv$obj@meta.data[[get_id_col(input$comp_ident%||%"cluster")]])
      }
      p <- p + guides(color=guide_legend(override.aes=list(size=4)))
      ggsave(f,fix_panel_size(p),width=dd$w,height=dd$h,dpi=600,bg="white",device="png")
    })

  output$dl_cu_ident <- downloadHandler(
    filename=function() paste0("comp_umap_identity_",Sys.Date(),".png"),
    content=function(f){
      obj <- apply_ident(rv$obj,input$comp_ident%||%"cluster")
      labs_v <- unique(as.character(Idents(obj)))
      dd  <- repel_dims(labs_v)
      obj <- filter_by_samples(obj, input$active_samples)
      p   <- umap_plot(obj,label=TRUE,pt.size=0.3,title="UMAP by Identity")+
             guides(color=guide_legend(override.aes=list(size=4),ncol=if(length(labs_v)>28) 2 else 1))+
             fixed_ident_scale(rv$obj@meta.data[[get_id_col(input$comp_ident%||%"cluster")]])
      ggsave(f,fix_panel_size(p),width=dd$w,height=dd$h,dpi=600,bg="white",device="png")
    })
  output$dl_comp_bar  <- downloadHandler(filename=function() paste0("composition_",Sys.Date(),".png"),content=function(f){id<-get_id_col(input$comp_ident);labs_v<-unique(as.character(rv$obj@meta.data[[id]]));dd<-barflip_dims(length(labs_v),labs_v);df<-rv$obj@meta.data%>%group_by(identity=.data[[id]],sample)%>%summarise(n=n(),.groups="drop")%>%group_by(identity)%>%mutate(pct=n/sum(n)*100);p<-ggplot(df,aes(reorder(identity,pct),pct,fill=sample))+geom_col(width=0.7)+coord_flip()+theme_minimal(base_size=13)+labs(x=NULL,y="%",fill="Sample")+theme(plot.margin=margin(t=10,r=10,b=10,l=max(10,max_label_len(labs_v)*1.5),unit="pt"));ggsave(f,p,width=dd$w,height=dd$h,dpi=600,bg="white")})
  output$dl_cg_feat   <- downloadHandler(filename=function() paste0("comp_feature_",Sys.Date(),".png"),content=function(f){req(input$comp_genes);gns<-resolve_genes(input$comp_genes,rv$obj);obj<-apply_ident(rv$obj,input$comp_ident%||%"cluster");obj<-filter_by_samples(obj,input$active_samples);n_samps<-length(unique(as.character(obj$sample)));n_genes<-length(gns);lims<-umap_full_xylim(rv$obj);p<-FeaturePlot(obj,features=gns,split.by="sample",reduction="umap",max.cutoff=3,pt.size=0.2,ncol=n_samps) & coord_cartesian(xlim=lims$x, ylim=lims$y);ggsave(f,p,width=max(6,5*n_samps),height=max(5,4*n_genes),dpi=600,bg="white")})
  output$dl_cg_vln <- downloadHandler(
    filename = function() {
      paste0(
        "comp_violin_",
        Sys.Date(),
        ".png"
      )
    },
    
    content = function(f) {
      req(input$comp_genes)
      
      gns <- resolve_genes(
        input$comp_genes,
        rv$obj
      )
      
      obj <- apply_ident(
        rv$obj,
        input$comp_ident %||% "cluster"
      )
      
      id_col <- get_id_col(input$comp_ident)
      
      clust_arg <- if(
        identical(input$cg_vln_scope %||% "global", "selected")
      ) {
        input$cg_vln_clusters_sel
      } else {
        input$cg_vln_clust
      }
      
      obj <- cg_scope_obj(
        obj,
        id_col,
        input$cg_vln_scope %||% "global",
        clust_arg
      )
      obj <- filter_by_samples(obj, input$active_samples)

      id_labels <- unique(
        as.character(obj@meta.data[[id_col]])
      )

      gns_use <- nonconstant_genes(obj, gns)
      req(length(gns_use) > 0)
      
      is_single <- identical(
        input$cg_vln_scope %||% "global",
        "percluster"
      )
      
      side_margin <- cg_vln_side_margin(id_labels)
      
      if (length(gns_use) == 1) {
        p <- VlnPlot(
          obj,
          features = gns_use,
          group.by = id_col,
          split.by = "sample",
          pt.size = 0
        ) +
          geom_point(
            position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.9),
            size = 0.3,
            alpha = 0.6
          ) &
          theme(
            plot.margin = margin(
              t = 60,
              r = side_margin,
              b = cg_vln_bottom_margin(id_labels),
              l = side_margin,
              unit = "pt"
            )
          )
      } else {
        p <- VlnPlot(
          obj,
          features = gns_use,
          group.by = id_col,
          split.by = "sample",
          stack = TRUE,
          flip = TRUE,
          pt.size = 0
        ) +
          geom_point(
            position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.9),
            size = 0.05,
            alpha = 0.15
          ) &
          theme(
            plot.margin = margin(
              t = 60,
              r = max(
                side_margin,
                max_label_len(gns_use) * 4
              ),
              b = cg_vln_bottom_margin(id_labels),
              l = side_margin,
              unit = "pt"
            ),
            panel.spacing = unit(0.6, "cm")
          )
      }
      
      p <- p &
        min_text_theme &
        scale_y_continuous(
          expand = expansion(mult = c(0.05, 0.15))
        ) &
        theme(
          axis.text.x = element_text(
            size = if (length(id_labels) <= 28) {
              13
            } else {
              pmax(
                9,
                13 - (length(id_labels) - 28) * 0.3
              )
            }
          ),
          axis.text.y = element_text(
            size = if (length(id_labels) <= 28) {
              13
            } else {
              pmax(
                9,
                13 - (length(id_labels) - 28) * 0.3
              )
            }
          ),
          axis.title = element_text(size = 12)
        )
      
      p <- cg_vln_apply_single(
        p,
        is_single,
        id_labels[1]
      )
      
      w_final <- if (is_single) {
        cg_vln_single_width_in
      } else {
        cg_vln_width_in(id_labels) * 1.12 + 1
      }
      
      h_final <- (
        max(6, length(gns_use) * 2.8) +
          cg_vln_bottom_margin(id_labels) / 72 +
          # When !is_single, ggplot's own rotated (45deg) identity tick
          # labels (Seurat's own VlnPlot default; cg_vln_apply_single only
          # resets this to horizontal when is_single) consume their own,
          # separate auto-sized vertical space in addition to the explicit
          # bottom plot.margin above -- without this second term, many long
          # identity labels make the two demands exceed the canvas height,
          # starving the panel's "null" row down to zero (confirmed via
          # ggplot_build()/gtable inspection). is_single mode never rotates
          # labels, so it never needs this extra budget.
          if (!is_single) cg_vln_bottom_margin(id_labels) / 72 else 0
      ) * 1.08 + 0.5

      # See cg_vln_width_in's comment: without this floor, is_single==FALSE
      # plots with a long identity label can starve the panel down to a
      # near-invisible sliver in both dimensions at once.
      if (!is_single) {
        h_final <- max(h_final, cg_vln_min_canvas_in(id_labels)$h)
      }

      ggsave(
        f,
        p,
        width = w_final,
        height = h_final,
        dpi = 600,
        bg = "white",
        limitsize = FALSE
      )
    }
  )
  output$dl_cg_dot    <- downloadHandler(filename=function() paste0("comp_dot_",Sys.Date(),".png"),content=function(f){req(input$comp_genes);gns<-resolve_genes(input$comp_genes,rv$obj);obj<-apply_ident(rv$obj,input$comp_ident%||%"cluster");id_col<-get_id_col(input$comp_ident);clust_arg<-if(identical(input$cg_dot_scope%||%"global","selected")) input$cg_dot_clusters_sel else input$cg_dot_clust;obj<-cg_scope_obj(obj,id_col,input$cg_dot_scope%||%"global",clust_arg);obj<-filter_by_samples(obj,input$active_samples);samp_names_use<-sort(unique(as.character(obj$sample)));ns<-length(samp_names_use);cols_u<-rep(c("blue","red","green","orange","purple"),ns)[1:ns];id_labels<-unique(as.character(obj@meta.data[[id_col]]));n_id<-length(id_labels);dd<-dotplot_dims(length(gns),n_id*ns,id_labels);p<-DotPlot(obj,features=gns,split.by="sample",cols=cols_u,dot.scale=8)+RotatedAxis()&theme(plot.margin=margin(t=10,r=10,b=max(60,max_label_len(id_labels)*3),l=10,unit="pt"))&min_text_theme&theme(axis.text.x=element_text(size=if(length(id_labels)<=28) 13 else pmax(9, 13-(length(id_labels)-28)*0.3)), axis.text.y=element_text(size=if(length(id_labels)<=28) 13 else pmax(9, 13-(length(id_labels)-28)*0.3)), axis.title=element_text(size=12));p<-with_dot_legend(p,samp_names_use,cols_u);ggsave(f,p,width=dd$w,height=dd$h,dpi=600,bg="white",device="png")})
  output$dl_pb        <- downloadHandler(filename=function() paste0("pseudobulk_",Sys.Date(),".png"),content=function(f){req(rv$pb_obj,input$pb_ct);s1<-get_s1();s2<-get_s2();ct<-as.character(input$pb_ct[1]);all_ids<-colnames(rv$pb_obj);i1<-find_pb_col(rv$pb_obj,ct,s1,rv$pb_id_col,rv$obj);i2<-find_pb_col(rv$pb_obj,ct,s2,rv$pb_id_col,rv$obj);if(is.na(i1)||is.na(i2)){png(f,width=800,height=800);plot.new();text(0.5,0.5,paste0("Column not found.\nAvailable:\n",paste(all_ids,collapse="\n")),cex=0.7);dev.off();return()};hi<-if(nchar(trimws(input$pb_hi))>0) resolve_genes(input$pb_hi,rv$obj) else character(0);r_val<-tryCatch({pb_mat<-GetAssayData(rv$pb_obj,assay=DefaultAssay(rv$pb_obj),layer="data");if(i1%in%colnames(pb_mat)&&i2%in%colnames(pb_mat)) round(cor(as.numeric(pb_mat[,i1]),as.numeric(pb_mat[,i2]),use="complete.obs"),3) else NA_real_},error=function(e) NA_real_);p<-CellScatter(rv$pb_obj,i1,i2,highlight=if(length(hi)) hi else NULL);if(length(hi)) p<-LabelPoints(plot=p,points=hi,repel=TRUE);p<-p+labs(title=paste0("Pseudobulk: ",ct,"  |  ",s1," vs ",s2),subtitle=if(!is.na(r_val)) paste0("Pearson r = ",r_val) else NULL)+theme(plot.subtitle=element_text(size=13,face="bold",color="#2c3e50",hjust=0.5));ggsave(f,p,width=8,height=8,dpi=600,bg="white")})
  output$dl_volc      <- dl600(function(){req(rv$de_res,rv$de_done);make_volcano(rv$de_res,lfc_thresh=as.numeric(input$volc_lfc)%||%1,pval_thresh=as.numeric(input$volc_pval)%||%0.05,top_n=as.integer(input$volc_n)%||%10,title=paste0("DE: ",input$de_clust%||%""," | log\u2082(",rv$de_s1," / ",rv$de_s2,")"),s1=rv$de_s1,s2=rv$de_s2)},"volcano",10,8)

  output$dl_gsea_plot <- downloadHandler(
    filename=function() paste0("gsea_",gsub("[^A-Za-z0-9]","_",rv$gsea_label%||%"plot"),"_",Sys.Date(),".png"),
    content=function(f){
      top_df<-gsea_top_df(); req(nrow(top_df)>0)
      dd<-barflip_dims(nrow(top_df),trunc_pw(top_df$pathway))
      dd$w <- dd$w + gsea_legend_extra_in  # external patchwork legend column, see gsea_with_legend()
      p<-gsea_make_plot(top_df,input$gsea_plot_type%||%"dot",rv$gsea_s1,rv$gsea_s2,rv$gsea_label)
      p<-gsea_with_legend(p)
      ggsave(f,p,width=dd$w,height=dd$h,dpi=600,bg="white",device="png")
    })
  output$dl_gsea_csv <- downloadHandler(
    filename=function() paste0("gsea_",gsub("[^A-Za-z0-9]","_",rv$gsea_label%||%"results"),"_",Sys.Date(),".csv"),
    content=function(f){
      req(rv$gsea_done,rv$gsea_res)
      df<-as.data.frame(rv$gsea_res)
      df$leadingEdge<-sapply(df$leadingEdge,function(x) paste(x,collapse="; "))
      write.csv(df,f,row.names=FALSE)
    })
  output$dl_pe_plot <- downloadHandler(
    filename=function() paste0("pathway_explorer_",gsub("[^A-Za-z0-9]","_",rv$gsea_label%||%"pe"),"_",Sys.Date(),".png"),
    content=function(f){
      req(rv$gsea_done,rv$gsea_res,input$pe_pathways,length(input$pe_pathways)>0)
      df<-rv$gsea_res%>%filter(pathway%in%input$pe_pathways)
      if(!nrow(df)) return()
      if(identical(input$pe_view%||%"barplot","trace")){
        req(!is.null(rv$gsea_pl),!is.null(rv$gsea_stats))
        plots<-lapply(input$pe_pathways,function(pw){
          if(pw%in%names(rv$gsea_pl))
            fgsea::plotEnrichment(rv$gsea_pl[[pw]],rv$gsea_stats)+
              labs(title=trunc_pw(pw,max_len=100))+theme_minimal(base_size=11)+min_text_theme
          else NULL
        })
        plots<-Filter(Negate(is.null),plots)
        if(!length(plots)) return()
        p<-wrap_plots(plots,ncol=1)
        ggsave(f,p,width=10,height=max(4,length(plots)*3.5),dpi=600,bg="white")
      } else {
        dd<-barflip_dims(nrow(df),trunc_pw(df$pathway))
        p<-pe_build_bar(df,rv$gsea_s1,rv$gsea_s2,rv$gsea_label)
        ggsave(f,p,width=dd$w,height=dd$h,dpi=600,bg="white")
      }
    })

  output$dl_ora_plot <- downloadHandler(
    filename=function() paste0("ora_",gsub("[^A-Za-z0-9]","_",rv$ora_label%||%"plot"),"_",Sys.Date(),".png"),
    content=function(f){
      df <- ora_combined()
      if(is.null(df)){showNotification("No significant ORA results to download.",type="warning");return()}
      dd <- barflip_dims(nrow(df), trunc_pw(df$Description))
      p  <- ora_make_plot(df, rv$ora_label, rv$ora_s1, rv$ora_s2, input$ora_plot_type%||%"bar")
      # device= is required here — Shiny hands downloadHandler an
      # extensionless temp path, which ggsave cannot infer a device from
      ggsave(f, p, width=dd$w, height=dd$h, dpi=600, bg="white", device="png", limitsize=FALSE)
    })

  output$dl_ora_csv <- downloadHandler(
    filename=function() paste0("ora_",gsub("[^A-Za-z0-9]","_",rv$ora_label%||%"results"),"_",Sys.Date(),".csv"),
    content=function(f){
      df <- ora_combined()
      if(is.null(df)){showNotification("No significant ORA results to download.",type="warning");return()}
      write.csv(df, f, row.names=FALSE)
    })

} # /server

shinyApp(ui=ui,server=server)

