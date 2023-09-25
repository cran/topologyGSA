.processGraph <- function(dag) {
  graph     <- graph_from_graphnel(dag)
  moral     <- moralize(graph)
  tg        <- triangulate(moral)
  adj.moral <- as(moral, "matrix")

  list(adj.moral = adj.moral,
       cli.moral = qpGetCliques(adj.moral, verbose=F),
       cli.tg    = getCliques(tg),
       moral     = moral,
       tg        = tg)
}
