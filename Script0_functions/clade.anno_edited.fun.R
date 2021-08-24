
clade.anno_edited <- function(gtree, anno.data, alpha=0.2, anno.depth=3, anno.x=10, anno.y=40){
  
  short.labs <- ""
  get_offset <- function(x) {(x*0.2+0.2)^2}
  get_angle <- function(node){
    data <- gtree$data
    sp <- tidytree::offspring(data, node)$node
    sp2 <- c(sp, node)
    sp.df <- data[match(sp2, data$node),]
    mean(range(sp.df$angle))
  }
  anno.data <- arrange(anno.data, node)
  hilight.color <- anno.data$color
  node_list <- anno.data$node
  node_ids <- (gtree$data %>% filter(label %in% node_list ) %>% arrange(label))$node
  anno <- rep('white', nrow(gtree$data))
  ## add hilight ... duplicated code
  ## FIXME: duplicated code...
  for(i in 1:length(node_ids)){
    n <- node_ids[i]
    color <- hilight.color[i]
    anno[n] <- color
    mapping <- gtree$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    offset <- get_offset(nodeClass)
    gtree <-
      gtree + geom_hilight(node=n, fill=color, alpha=alpha,
                           extend=offset)
  }
  gtree$layers <- rev(gtree$layers)
  gtree <- gtree + geom_point2(aes(size=I(nodeSize)), fill=anno, shape=21)
  ## add labels
  short.labs.anno <- NULL
  for(i in 1:length(node_ids)){
    n <- node_ids[i]
    mapping <- gtree$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    if(nodeClass <= anno.depth){## species and strains
      lab <- short.labs[1]
      short.labs <- short.labs[-1]
      # short.labs.anno <- paste0(short.labs.anno, sep='\n', paste0(lab, ': ', mapping$label))
      if(is.null(short.labs.anno)){
        short.labs.anno = data.frame(lab=lab, annot = mapping$label, stringsAsFactors = F)
      }else{
        short.labs.anno = rbind(short.labs.anno,
                                c(lab, mapping$label))
      }
    }
    else{
      lab <- mapping$label
    }
    offset <- get_offset(nodeClass) - 0.4
    angle <- get_angle(n) + 90
    gtree <- gtree +
      geom_cladelabel(node=n, label=lab, angle=angle,
                      fontsize=1.5+sqrt(nodeClass),
                      offset=offset, barsize=0, hjust=0.5)
  }
  if(is.null(short.labs.anno)){return(gtree)}
  ## add short labels
  anno_shapes = sapply(short.labs.anno$lab, utf8ToInt)
  
  #(added by Gabriel 210823) ####
  set_hilight_legend <- function(p, color, label, alpha=.4) {
    d <- data.frame(color=factor(color,labels = c("green","orange", "indianred3", 'darkgreen', "darkblue", "darkred", "gray80")), 
                    var=factor(label,labels=c("AHI", "T90", "BMI", "AHI+T90", "AHI+BMI", "T90+BMI", "All")), 
                    x=0, y=1, alpha=alpha)
    p + geom_rect(aes_(fill=~var, xmin=~x, xmax=~x, ymin=~y, ymax=~y), data=d, inherit.aes=F) +
      guides(fill=guide_legend(override.aes=list(fill=alpha(c("green","orange", "indianred3", 'darkgreen', "darkblue", "darkred", "gray80"),
      d$alpha)), title=NULL))
      #guides(fill=guide_legend(override.aes=list(fill=alpha(d$color, d$alpha)), title=NULL))
  }
  
  cls = c("green","orange", "indianred3", 'darkgreen', "darkblue", "darkred", "gray80")
  names = c("AHI", "T90", "BMI", "AHI+T90", "AHI+BMI", "T90+BMI", "All")
  set_hilight_legend(p= gtree, color = cls, label = names) + theme(legend.position="right")
  
  
  #### 
  

  
  
  #gtree + geom_point(data = short.labs.anno,
   #                  aes(x=0, y=0, shape = factor(annot)),
    #                 size=0, stroke=0) +
#    guides(
#      shape = guide_legend(override.aes = list(size=3, shape=anno_shapes))
#    ) +
#    theme(legend.position = 'none',#c(1.2,0.5),
 #         legend.title = element_blank())
}