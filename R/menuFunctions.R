
## set out menu logic for downstream

.menuItemLogic <- function(flags) {
	## initiate the menue logic for PASS, WARN, FAIL flags
	## menuLogic will be a list of length 5 containing information on flags
	menuLogic <- list()
	if (all(flags$Status == "PASS")) {
		menuLogic[[1]] <- "PASS"
		menuLogic[[2]] <- "green"
	}
	else {
		menuLogic[[1]] <- "WARN"
		menuLogic[[2]] <- "yellow"
	}
	
	if (all(flags$Status == "FAIL")) {
		menuLogic[[1]] <- "FAIL"
		menuLogic[[2]] <- "red"
	}
	
	# Fail values
	menuLogic[[3]] <- c(sum(flags$Status == "FAIL"), length(flags$Status))
	
	menuLogic[[4]] <- c(sum(flags$Status == "WARN"), length(flags$Status))
	
	# Pass values
	menuLogic[[5]] <-c(sum(flags$Status == "PASS"), length(flags$Status))
	
	menuLogic
}

  ###### renders the box for presence or access of an icon
.renderValBox <- function(count, status, ic, c) {
	renderValueBox({
		valueBox(
			value = paste(count[1], count[2], sep = "/"),
			subtitle = status,
			icon = icon(ic, class = "fa-lg"),
			color = c
		)
	})
}