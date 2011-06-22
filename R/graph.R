.processGraph <- function(dag) {
  moral     <- moralize(dag)
  tg        <- triangulate(moral)
  adj.moral <- as(moral, "matrix")

  list(adj.moral = adj.moral,
       cli.moral = qpGetCliques(adj.moral, verbose=F),
       cli.tg    = maxClique(tg),
       moral     = moral,
       tg        = tg)
}
