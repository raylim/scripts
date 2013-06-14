#!/usr/bin/env Rscript
# runs rsp on a file, passes on arguments

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library("R.rsp"));

allArguments <- commandArgs(trailingOnly = T)

if (length(allArguments) < 1) {
    cat("Need rsp file");
    stop();
}

rspFile <- allArguments[1];
arguments <<- allArguments[-1];


rsp.default <- structure(function (filename = NULL, path = NULL, text = NULL, 
    response = NULL, ..., envir = parent.frame(), postprocess = TRUE, 
    verbose = FALSE) 
{
    suppressPackageStartupMessages(require("R.rsp", quietly = TRUE)) || 
        throw("Package not loaded: R.rsp")
    rVer <- as.character(getRversion())
    if (compareVersion(rVer, "2.13.0") < 0) {
        tempfile <- function(..., fileext = "") {
            for (kk in 1:100) {
                pathnameT <- base::tempfile(...)
                if (fileext != "") {
                  pathnameT <- sprintf("%s%s", pathnameT, fileext)
                  if (!file.exists(pathnameT)) {
                    return(pathnameT)
                  }
                }
            }
            stop("Failed to create a non-existing temporary pathname.")
        }
    }
    rspPlain <- function(pathname, response = NULL, ..., verbose = FALSE) {
        if (is.null(response)) {
            verbose && enter(verbose, "Creating FileRspResponse")
            pattern <- "((.*)[.]([^.]+))[.]([^.]+)$"
            pathname2 <- gsub(pattern, "\\1", pathname)
            pathname2 <- Arguments$getWritablePathname(pathname2)
            response <- FileRspResponse(file = pathname2, overwrite = TRUE)
            verbose && exit(verbose)
        }
        if (inherits(response, "connection")) {
            response <- FileRspResponse(file = response)
        }
        else if (is.character(response)) {
            pathname <- Arguments$getWritablePathname(response)
            response <- FileRspResponse(file = pathname)
        }
        else if (!inherits(response, "RspResponse")) {
            throw("Argument 'response' is not an RspResponse object: ", 
                class(response)[1])
        }
        verbose <- Arguments$getVerbose(verbose)
        if (verbose) {
            pushState(verbose)
            on.exit(popState(verbose))
        }
        verbose && enter(verbose, "Compiling RSP-embedded plain document")
        verbose && cat(verbose, "Input pathname: ", pathname)
        verbose && printf(verbose, "%s:\n", class(response)[1])
        verbose && print(verbose, response)
        pathname2 <- getOutput(response)
        verbose && cat(verbose, "Response output class: ", class(pathname2)[1])
        verbose && cat(verbose, "Response output pathname: ", 
            pathname2)
        verbose && enter(verbose, "Calling sourceRspV2()")
        sourceRspV2(pathname, path = NULL, ..., response = response, 
            envir = envir, verbose = less(verbose, 20))
        verbose && exit(verbose)
        verbose && exit(verbose)
        invisible(pathname2)
    }
    if (!is.null(filename) && !is.null(text)) {
        throw("Only one of arguments 'filename' and 'text' can be specified.")
    }
    if (!is.null(text)) {
        pathnameT <- tempfile(fileext = ".txt.rsp")
        on.exit(file.remove(pathnameT))
        writeLines(text = text, con = pathnameT)
        filename <- pathnameT
        path <- NULL
        if (is.null(response)) {
            response <- stdout()
        }
    }
    pathname <- Arguments$getReadablePathname(filename, path = path, 
        mustExist = TRUE)
    postprocess <- Arguments$getLogical(postprocess)
    verbose <- Arguments$getVerbose(verbose)
    if (verbose) {
        pushState(verbose)
        on.exit(popState(verbose))
    }
    verbose && enter(verbose, "Compiling RSP document")
    verbose && cat(verbose, "RSP pathname: ", pathname)
    pattern <- "((.*)[.]([^.]+))[.]([^.]+)$"
    ext <- gsub(pattern, "\\3", pathname)
    type <- tolower(ext)
    verbose && cat(verbose, "RSP type: ", type)
    verbose && cat(verbose, "Postprocess (if recognized): ", 
        postprocess)
    postProcessor <- NULL
    if (postprocess) {
        verbose && enter(verbose, "Searching for document-type specific postprocessor")
        postProcessors <- list(tex = compileLaTeX, rnw = compileSweave)
        for (key in names(postProcessors)) {
            pattern <- key
            if (regexpr(pattern, type) != -1) {
                postProcessor <- postProcessors[[key]]
                verbose && cat(verbose, "Match: ", key)
                break
            }
        }
        if (is.null(postProcessor)) {
            verbose && cat(verbose, "Postprocessor found: <none>")
        }
        else {
            verbose && cat(verbose, "Postprocessor found: ", 
                type)
        }
        verbose && exit(verbose)
    }
    verbose && enter(verbose, "Preprocessing, translating, and evaluating RSP document")
    res <- rspPlain(pathname, response = response, ..., verbose = verbose)
    wasFileGenerated <- inherits(res, "character")
    if (wasFileGenerated) {
        pathname2 <- res
        verbose && cat(verbose, "Output pathname: ", pathname2)
    }
    verbose && exit(verbose)
    if (!is.null(postProcessor)) {
        if (wasFileGenerated) {
            verbose && enter(verbose, "Postprocessing generated document")
            verbose && cat(verbose, "Input pathname: ", pathname2)
            pathname3 <- postProcessor(pathname2, ..., verbose = verbose)
            verbose && cat(verbose, "Output pathname: ", pathname3)
            verbose && exit(verbose)
            res <- pathname3
        }
    }
    if (wasFileGenerated) {
        verbose && cat(verbose, "Output document pathname: ", 
            res)
    }
    else {
        verbose && printf(verbose, "Output written to: %s [%d]\n", 
            class(res)[1], res)
    }
    verbose && exit(verbose)
    invisible(res)
}, modifiers = "public")

rsp(rspFile, verbose = T, postprocess = F);
