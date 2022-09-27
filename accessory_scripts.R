seqToHash64Base62 <- function(seqs) {
  sapply(
    seqs,
    function(x) Rmpfr::formatMpfr(
      Rmpfr::mpfr(
        digest::digest(x, "xxhash64"),
        base = 16,
        precBits = 64
      ),
      base = 62,
      format = "d",
      decimal.mark = "",
      digits = 0,
      drop0trailing = F,
      mode = "integer"
    )
  )
}
