loadBAMMtools <- function(){
	scripts <- list.files()
	scripts <- grep(pattern = '.R',scripts, value=T, fixed=T)
	for (i in 1:length(scripts)){
		source(scripts[i])
	}
}