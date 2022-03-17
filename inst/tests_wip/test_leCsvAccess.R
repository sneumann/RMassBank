# Test correct results of findMz.formula wiht positive and negative charge,
# single and multiple charge, no charge, and fictitious negative atoms
expect_equal(findMz.formula("C6", "")$mzCenter, 72)
expect_equal(
  findMz.formula("C6", "mH")$mzCenter, 
  72 - 1.0078 + RMassBank:::.emass,
  tolerance = 0.00001 )
expect_equal(
  findMz.formula("C6", "pH")$mzCenter, 
  72 + 1.0078 - RMassBank:::.emass,
  tolerance = 0.00001 )
expect_equal(
  findMz.formula("C6H-1", "")$mzCenter, 
  72 - 1.0078,
  tolerance = 0.00001 )
expect_equal(
  findMz.formula("C6H-1", "mM")$mzCenter, 
  72 - 1.0078 + RMassBank:::.emass,
  tolerance = 0.00001 )
expect_equal(
  findMz.formula("C6H-1", "pM")$mzCenter, 
  72 - 1.0078 - RMassBank:::.emass,
  tolerance = 0.00001 )
expect_equal(
  findMz.formula("C6", "m2H_c2")$mzCenter, 
  (72 - (2*1.0078) + 2*RMassBank:::.emass) / 2,
  tolerance = 0.00001 )
expect_equal(
  findMz.formula("C6H-1", "m2H_c2")$mzCenter, 
  (72 - (3*1.0078) + 2*RMassBank:::.emass) / 2,
  tolerance = 0.00001 )




