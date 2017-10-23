## File Name: rm_sdt_probraterfct1.R
## File Version: 0.02
## File Last Change: 2017-10-02 17:45:54


################################################################
# C code
# calculation of probabilities
rm_sdt_probraterfct1 <- function(crater,drater,dimA,B,dimB){ 
	rm_probraterfct1(CRA=crater, DRA=drater, dimAA=dimA, BB=B, dimBB=dimB)
}

probraterfct1 <- rm_sdt_probraterfct1
