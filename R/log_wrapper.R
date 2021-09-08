#' @import logger
NULL

#' Pass arguments to logger::log_info using custom RMassBank-logging settings
#'
#' The logging file to be used can be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_info
#' @author pstahlhofen
#' @export
rmb_log_info <- function(...) {
	logger::log_info(..., namespace='RMassBank')
}

#' Pass arguments to logger::log_trace using custom RMassBank-logging settings
#'
#' The logging file to be used can be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_trace
#' @author pstahlhofen
#' @export
rmb_log_trace <- function(...) {
	logger::log_trace(..., namespace='RMassBank')
}

#' Pass arguments to logger::log_debug using custom RMassBank-logging settings
#'
#' The logging file to be used can be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_debug
#' @author pstahlhofen
#' @export
rmb_log_debug <- function(...) {
	logger::log_debug(..., namespace='RMassBank')
}

#' Pass arguments to logger::log_warn using custom RMassBank-logging settings
#'
#' The logging file to be used can be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_warn
#' @author pstahlhofen
#' @export
rmb_log_warn <- function(...) {
	logger::log_warn(..., namespace='RMassBank')
}

#' Pass arguments to logger::log_success using custom RMassBank-logging settings
#'
#' The logging file to be used can be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_success
#' @author pstahlhofen
#' @export
rmb_log_success <- function(...) {
	logger::log_success(..., namespace='RMassBank')
}

#' Pass arguments to logger::log_error using custom RMassBank-logging settings
#'
#' The logging file to be used can be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_error
#' @author pstahlhofen
#' @export
rmb_log_error <- function(...) {
	logger::log_error(..., namespace='RMassBank')
}

#' Pass arguments to logger::log_fatal using custom RMassBank-logging settings
#'
#' The logging file to be used can be specified by the user in the \code{logging_file} field of \code{settings.ini}
#' @seealso logger::log_fatal
#' @author pstahlhofen
#' @export
rmb_log_fatal <- function(...) {
	logger::log_fatal(..., namespace='RMassBank')
}
