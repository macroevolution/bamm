assignColorBreaks = function(NCOLORS = 64, spex = "s") {
	if (!exists("rates", envir = .dtRatesEnv)) {
		stop("Could not find 'rates' in .dtRatesEnv");
	}
	dtr = get("rates", envir = .dtRatesEnv);
	if (mode(dtr$rates) == "numeric") {
		bks = quantile(dtr$rates, seq(0,1, length.out=(NCOLORS+1)));
		assign("colorbreaks", bks, envir = .dtRatesEnv);
	}
	else if(mode(dtr$rates) == "list") {
		if (tolower(spex) == "s") {
			bks = quantile(dtr$rates[[1]], seq(0,1, length.out=(NCOLORS+1)));
		}
		else if (tolower(spex) == "e") {
			bks = quantile(dtr$rates[[2]], seq(0,1, length.out=(NCOLORS+1)));
		}
		else {
			bks = quantile(dtr$rates[[1]] - dtr$rates[[2]], seq(0,1, length.out=(NCOLORS+1)));
		}
		assign("colorbreaks", bks, envir = .dtRatesEnv);
		assign("colorscall", get("call", .dtRatesEnv), envir = .dtRatesEnv);
	}
}