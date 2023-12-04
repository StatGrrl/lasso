color.tint <- function(color, alpha=255, tint=0) {
  # INPUT:
  # color...color name, see color()
  # alpha...tranparency
  #         integer between 0 (complete) and 255 (none, default)
  # tint....how much to lighten (< 0) or darken (> 0)
  #         float between -1 (white) and 1 (black)
  #         default is 0 (none)
  # OUTPUT:
  # returns tinted color
  if (alpha < 0 | alpha > 255) stop("alpha must be between 0 and 255")
  if (abs(tint) > 1) stop("t must be between -1 and 1")
  col <- col2rgb(color)
  if (tint > 0) new <- col * (1 - tint) # darken
  if (tint < 0) new <-  col + ((-tint) * (255 - col)) # lighten
  if (tint==0) {
    out <- rgb(col[1], col[2], col[3], alpha, max=255)
  } else out <- rgb(new[1], new[2], new[3], alpha, max=255)
  out
}

color.vec <- function(vec, color, alpha=255) {
  # INPUT:
  # vec.....vector of proportions
  # color...color name, see color()
  # alpha...tranparency
  #         integer between 0 (complete) and 255 (none, default)
  #
  # OUTPUT:
  # vector of colors for vec
  # legend or key data frame for colors
  
  key <- data.frame(from=seq(0, 0.95, by=0.05), 
                     to=seq(0.05, 1, by=0.05), 
                     tint=seq(-1, 0.9, by=0.1),
                     label=rep(0,20), color=rep(0,20))
  for (i in 1:nrow(key)) {
    key[i,"label"] <- paste(key[i,"to"]*100, "%", sep="")
    key[i,"color"] <- color.tint(color, alpha, key[i,"tint"])
  }
  vec.pos <- unlist(lapply(vec, findInterval, key[,"from"]))
  vec.color <- key[vec.pos,"color"]
  list(vec=vec.color, key=key)
}