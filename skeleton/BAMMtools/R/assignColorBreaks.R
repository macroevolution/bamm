assignColorBreaks <- function(rates, NCOLORS = 64, spex = "s") {
	if (mode(rates) == "numeric") {
		bks <- quantile(rates, seq(0,1, length.out=(NCOLORS+1)));
	}
	else if(mode(rates) == "list") {
		if (tolower(spex) == "s") {
			bks <- quantile(rates[[1]], seq(0,1, length.out=(NCOLORS+1)));
		}
		else if (tolower(spex) == "e") {
			bks <- quantile(rates[[2]], seq(0,1, length.out=(NCOLORS+1)));
		}
		else {
			bks <- quantile(rates[[1]] - rates[[2]], seq(0,1, length.out=(NCOLORS+1)));
		}
	}
	return(bks);
}
