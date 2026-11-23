############################################################
## Figure S8
##
## a) MAG-weighted observed categories versus the genus null
## b) Observed categories under MAG and species weighting
## c) Exact paired classification transitions
## d) Mixed fractions only: observed versus the two matched nulls
##
############################################################

library(grid)

figure_width <- 12
figure_height <- 9
# png_dpi <- 600

dir.create("figures", recursive = TRUE, showWarnings = FALSE)
dir.create("../envs", recursive = TRUE, showWarnings = FALSE)

status_levels <- c("Conserved present", "Mixed", "Conserved absent")
weightings <- c("MAG weighted", "Equal species")
status_colors <- c("Conserved present" = "#2166AC", "Mixed" = "#F4A582", "Conserved absent" = "#D9D9D9")
ink <- "black"
grid_color <- "#E5E5E5"

counts <- read.delim("sources/FigS8_status_counts.tsv",check.names=FALSE)
transitions <- read.delim("sources/FigS8_transitions.tsv",check.names=FALSE)
nulls <- read.delim("sources/FigS8_null_summary.tsv",check.names=FALSE)
paired <- read.delim("sources/FigS8_paired_classifications.tsv",check.names=FALSE)

N <- nrow(paired)
unchanged <- sum(paired$MAG_weighted_status == paired$equal_species_status)

counts$percent <- counts$n_combinations/N*100
M <- matrix(0,3,3,dimnames=list(status_levels,status_levels))
for(i in seq_len(nrow(transitions))) M[transitions$MAG_weighted_status[i],transitions$equal_species_status[i]] <- transitions$n_combinations[i]

row_pct <- 100*sweep(M,1,rowSums(M),"/")

# font size
axis_text_size <- 10.5
axis_title_size <- 12
legend_text_size <- 10
tag_size <- 14

txt <- function(label,x,y,size=axis_text_size,color=ink,font=1,just="left",rot=0) {
  grid.text(label,x=x,y=y,just=just,rot=rot,
    gp=gpar(fontsize=size,col=color,fontface=font,fontfamily="sans",lineheight=1.05))
}
seg <- function(x,y,col=ink,lwd=0.6,lty=1) grid.lines(x,y,gp=gpar(col=col,lwd=lwd,lty=lty))
rect <- function(x,y,w,h,fill="white",border=ink,lwd=0.6) grid.rect(x,y,w,h,gp=gpar(fill=fill,col=border,lwd=lwd))
pt <- function(x,y,pch=21,fill=ink,size=0.102) grid.points(x,y,pch=pch,size=unit(size,"in"),gp=gpar(fill=fill,col=ink,lwd=0.8))
tag <- function(s) txt(s,0.10,0.9,size=tag_size,font=1)

plot_frame <- function(x0,y0,w,h,xmax,ticks,rows=numeric(0)) {
  rect(x0+w/2,y0+h/2,w,h,fill="white",border=NA)
  for(t in ticks) seg(rep(x0+w*t/xmax,2),c(y0,y0+h),col=grid_color,lwd=0.35)
  for(y in rows) seg(c(x0,x0+w),c(y,y),col=grid_color,lwd=0.35)
  rect(x0+w/2,y0+h/2,w,h,fill=NA,border=ink,lwd=0.6)
  for(t in ticks) {
    x <- x0+w*t/xmax
    seg(c(x,x),c(y0,y0-0.013))
    txt(paste0(t,"%"),x,y0-0.049,size=axis_text_size,just="centre")
  }
}
comparison_legend <- function(y=0.80) {
  pt(0.32,y,pch=21,fill=ink,size=0.09)
  txt("Observed",0.344,y,size=legend_text_size)
  pt(0.60,y,pch=23,fill="white",size=0.095)
  txt("Null mean",0.624,y,size=legend_text_size)
}
null_comparison <- function(y,observed,n,fill,x0,width,xmax,label_p=FALSE) {
  xp <- function(z) x0+width*z/xmax
  seg(c(xp(observed),xp(n$null_mean_pct)),c(y,y),col="#ACACAC",lwd=0.65)
  seg(c(xp(n[["null_2.5_pct"]]),xp(n[["null_97.5_pct"]])),c(y,y),lwd=1)
  for(e in c(n[["null_2.5_pct"]],n[["null_97.5_pct"]])) seg(rep(xp(e),2),c(y-0.022,y+0.022),lwd=0.6)
  pt(xp(observed),y,pch=21,fill=fill,size=0.105)
  pt(xp(n$null_mean_pct),y,pch=23,fill="white",size=0.105)
  txt(sprintf("%.1f",observed),xp(observed),y+0.06,size=10.5,just="centre")
  txt(sprintf("%.1f",n$null_mean_pct),xp(n$null_mean_pct),y+0.06,size=10.5,just="centre")
  if(label_p) txt("P < 0.001",(xp(observed)+xp(n$null_mean_pct))/2,y-0.061,size=10,font=3,just="centre")
}

panel_a <- function() {
  tag("a")
  comparison_legend()
  x0 <- 0.30; y0 <- 0.20; width <- 0.66; height <- 0.54
  ys <- c(0.65,0.47,0.29)
  plot_frame(x0,y0,width,height,80,seq(0,80,20),rows=ys)
  for(i in 1:3) {
    s <- status_levels[i]
    o <- counts$percent[counts$weighting=="MAG weighted" & counts$status==s]
    n <- nulls[nulls$weighting=="MAG weighted" & nulls$status==s,]
    txt(c("Conserved\npresent","Mixed","Conserved\nabsent")[i],x0-0.024,ys[i],just="right")
    null_comparison(ys[i],o,n,status_colors[s],x0,width,80,label_p=s=="Mixed")
  }
  txt("Genus-pathway combinations (%)",x0+width/2,0.073,size=axis_title_size,just="centre")
}

panel_b <- function() {
  tag("b")
  x0 <- 0.27; width <- 0.69; ys <- c(0.60,0.36)
  plot_frame(x0,0.20,width,0.54,100,seq(0,100,20),rows=ys)
  for(i in 1:2) {
    w <- weightings[i]
    txt(c("MAG-weighted","Species-weighted")[i],x0-0.027,ys[i],just="right")
    cursor <- x0
    for(j in 1:3) {
      s <- status_levels[j]
      value <- counts$percent[counts$weighting==w & counts$status==s]
      delta <- width*value/100
      rect(cursor+delta/2,ys[i],delta,0.12,fill=status_colors[s],border="white",lwd=0.6)
      txt(sprintf("%.1f",value),cursor+delta/2,ys[i],size=10.5,color=if(j==1) "white" else ink,just="centre")
      cursor <- cursor+delta
    }
  }
  txt("Genus-pathway combinations (%)",x0+width/2,0.073,size=axis_title_size,just="centre")
}

panel_c <- function() {
  tag("c")
  x0 <- 0.33; y0 <- 0.20; cell_w <- 0.21; cell_h <- 0.18
  short <- c("Conserved\npresent","Mixed","Conserved\nabsent")
  cols <- colorRampPalette(c("#F7FBFF","#2166AC"))(101)
  for(j in 1:3) txt(short[j],x0+(j-0.5)*cell_w,0.137,size=10.5,just="centre")
  for(i in 1:3) {
    y <- y0+(3-i+0.5)*cell_h
    txt(short[i],x0-0.025,y,size=10.5,just="right")
    for(j in 1:3) {
      x <- x0+(j-0.5)*cell_w
      value <- row_pct[i,j]
      rect(x,y,cell_w,cell_h,fill=cols[round(value)+1],border="white",lwd=1.4)
      txt(paste0(format(M[i,j],big.mark=","),"\n",sprintf("%.1f%%",value)),x,y,size=11,
          color=if(value>60) "white" else ink,just="centre")
    }
  }
  rect(x0+1.5*cell_w,y0+1.5*cell_h,3*cell_w,3*cell_h,fill=NA,lwd=0.6)
  txt("MAG-weighted classification",0.045,0.47,size=axis_title_size,just="centre",rot=90)
  txt("Species-weighted classification",x0+1.5*cell_w,0.045,size=axis_title_size,just="centre")
}

panel_d <- function() {
  tag("d")
  comparison_legend()
  x0 <- 0.27; width <- 0.69; ys <- c(0.60,0.36)
  plot_frame(x0,0.20,width,0.54,50,seq(0,50,10),rows=ys)
  for(i in 1:2) {
    w <- weightings[i]
    txt(c("MAG-weighted","Species-weighted")[i],x0-0.027,ys[i],just="right")
    o <- counts$percent[counts$weighting==w & counts$status=="Mixed"]
    n <- nulls[nulls$weighting==w & nulls$status=="Mixed",]
    null_comparison(ys[i],o,n,status_colors["Mixed"],x0,width,50,label_p=TRUE)
  }
  txt("Mixed genus-pathway combinations (%)",x0+width/2,0.073,size=axis_title_size,just="centre")
}

draw_figure <- function() {
  grid.newpage()
  txt("Trait status",0.245,0.965,size=11)
  for(j in 1:3) {
    x <- c(0.345,0.55,0.66)[j]
    rect(x,0.965,0.017,0.022,fill=status_colors[j],border=NA)
    txt(status_levels[j],x+0.013,0.965,size=legend_text_size)
  }
  pushViewport(viewport(x=0.50,y=0.48,width=0.95,height=0.88,layout=grid.layout(2,2)))
  panels <- list(panel_a,panel_b,panel_c,panel_d)
  for(i in 1:4) {
    pushViewport(viewport(layout.pos.row=if(i<=2) 1 else 2,layout.pos.col=if(i%%2==1) 1 else 2))
    panels[[i]]()
    popViewport()
  }
  popViewport()
}

# Capture the complete figure
p_combined <- grid::grid.grabExpr(
  draw_figure(),
  width = figure_width,
  height = figure_height
)

# Display in the RStudio plot window
grid::grid.newpage()
grid::grid.draw(p_combined)

# Save PDF
ggplot2::ggsave(
  "figures/FigS8.pdf",
  plot = p_combined,
  width = figure_width,
  height = figure_height,
  units = "in"
)

# Optional PNG export
# ggplot2::ggsave(
#   "figures/FigS8.png",
#   plot = p_combined,
#   width = figure_width,
#   height = figure_height,
#   units = "in",
#   dpi = 600
# )

writeLines(
  capture.output(sessionInfo()),
  "../envs/FigS8.sessionInfo.txt"
)
