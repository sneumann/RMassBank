\name{smiles2mass}
\alias{smiles2mass}
\title{Calculates the mass from the smiles formula}
\usage{
  smiles2mass(SMILES)
}
\arguments{
  \item{SMILES}{The SMILES-string}
}
\value{
  A numeric with the calculated mass
}
\description{
  This function calculates the exact mass from a SMILES-String and is needed for a test
}
\examples{
smiles2mass(c(ccc1C(=NN=C2C)C(=O)N2)cc1)

}
\author{
  Erik Mueller
}

