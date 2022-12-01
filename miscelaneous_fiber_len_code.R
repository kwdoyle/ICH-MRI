# for (datype in datnames) {
#   # plot names
#   plt_names <- paste(datype, "_", consc_groups, sep="")
#   # assign plots
#   # the length of plot name variable SHOULD match the length of consc_groups
#   # because, y'know, they were used to make them in the first place..
#   # so I can just iterate over plt_names and use that index to pull out
#   # the corresponding consc_group too.
#   for (i in 1:length(plt_names)) {
#     assign(
#       plt_names[i],
#       ggplot()
#     )
#   }
# }



# for each data type table, make a plot for each conscious group and save it to a variable
for (name in names(parsedData)) {
  # make plot names
  plt_names <- paste(name, "_", consc_groups, sep="")
  # assign plots
  # the length of plot name variable SHOULD match the length of consc_groups
  # because, y'know, they were used to make them in the first place..
  # so I can just iterate over plt_names and use that index to pull out
  # the corresponding consc_group too.
  for (i in 1:length(plt_names)) {
    assign(
      plt_names[i],
      # need to use aes_string instead to pass each consc_group value each time.
      ggplot(parsedData[[name]], aes_string(x="len", y="mean", group="ID", color=consc_groups[i])) +
        geom_line(size=1.1) + theme_bw() +
        labs(title=paste(name, consc_groups[i]), x="length", y="mean", color="key")
    )
  }
}





plot_groups <- c("conscious_throughout", "recove")

PlotCombinedGroups <- function(datalist, control=FALSE) {
  out <- list()
  names_to_use <- names(datalist)[names(datalist) != "MainTable"]
  
  
  
  
}