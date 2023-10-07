#!/opt/conda/bin/Rscript

library("R6")
library("jsonlite")

#' 创建命令行，参照argparse
#' 
#' @export
ArgumentParser <- function(..., python_cmd = NULL,Parsername = "parser") { # nolint
  python_cmd <- find_python_cmd(python_cmd)
  assert_python_cmd(python_cmd)
  initial_python_code <- c("import argparse",
                           "try:",
                           "    import json",
                           "except ImportError:",
                           "    import simplejson as json",
                           "",
                           "def logical(s):",
                           "    if isinstance(s, bool):",
                           "        return bool",
                           "    elif s in ('T', 'TRUE', 'True', 'true'):",
                           "        return True",
                           "    elif s in ('F', 'FALSE', 'False', 'false'):",
                           "        return False",
                           "    else:",
                           "        raise ValueError(\"could not convert string to logical: '{}'\".format(s))",
                           "",
                           sprintf("%s = argparse.ArgumentParser(%s)",
                                   Parsername,
                                   convert_..._to_arguments("ArgumentParser", ...)),
                           "")
  Parser$new(Code$new(python_cmd, initial_python_code), Parsername)
}

Code <- R6Class("Code", # nolint
                public = list(
                  append = function(new_code) {
                    private$code <- c(private$code, new_code)
                  },
                  run = function(new_code) {
                    code <- c(private$code, new_code)
                    suppressWarnings(system2(private$cmd, input = code,
                                             stdout = TRUE, stderr = TRUE))
                  },
                  initialize = function(cmd, code = character(0)) {
                    private$cmd <- cmd
                    private$code <- code
                  }
                ),
                private = list(code = NULL, cmd = NULL)
)

Group <- R6Class("Group", # nolint
                 public = list(
                   initialize = function(python_code, name) {
                     private$python_code <- python_code
                     private$name <- name
                   },
                   add_argument = function(...,Parsername = private$name) {
                     private$python_code$append(sprintf("%s.add_argument(%s)",
                                                        Parsername, convert_..._to_arguments("add_argument", ...)))
                     return(invisible(NULL))
                   }
                 ),
                 private = list(python_code = NULL, name = NULL)
)

Subparsers <- R6Class("Subparsers", # nolint
                      public = list(
                        initialize = function(python_code, name) {
                          private$python_code <- python_code
                          private$name <- name
                        },
                        add_parser = function(...,Parsername = private$name) {
                          parser_name <- paste0(Parsername, "_subparser", private$n_subparsers)
                          private$n_subparsers <- private$n_subparsers + 1
                          private$python_code$append(sprintf("%s = %s.add_parser(%s)",
                                                             parser_name, Parsername,
                                                             convert_..._to_arguments("ArgumentParser", ...)))
                          Parser$new(private$python_code, parser_name)
                        }
                      ),
                      private = list(python_code = NULL, name = NULL, n_subparsers = 0)
)

Parser <- R6Class("Parser", # nolint
                  public = list(
                    parse_args = function(args = commandArgs(TRUE)) {
                      new_code <- c(sprintf("args = %s.parse_args([%s])", private$name,
                                            paste(sprintf("'%s'", args), collapse = ", ")),
                                    "print(json.dumps(args.__dict__, sort_keys=True))")
                      output <- private$python_code$run(new_code)
                      parse_args_output(output)
                    },
                    parse_known_args = function(args = commandArgs(TRUE)) {
                      new_code <- c(sprintf("args_remainder = %s.parse_known_args([%s])", private$name,
                                            paste(sprintf("'%s'", args), collapse = ", ")),
                                    "print(json.dumps((args_remainder[0].__dict__, args_remainder[1])))")
                      output <- private$python_code$run(new_code)
                      parse_args_output(output)
                    },
                    print_help = function() {
                      cat(private$python_code$run(sprintf("%s.print_help()", private$name)), sep = "\n")
                      invisible(NULL)
                    },
                    print_usage = function() {
                      cat(private$python_code$run(sprintf("%s.print_usage()", private$name)), sep = "\n")
                      invisible(NULL)
                    },
                    newparser = function(...,Parsername = paste0(private$name, private$n_parsers)) {
                      private$python_code$append(sprintf("%s = argparse.ArgumentParser(%s)",
                                                         Parsername,
                                                         convert_..._to_arguments("ArgumentParser", ...)))
                      private$n_parsers <- private$n_parsers+1
                      Parser$new(private$python_code, 
                                 Parsername,
                                 private$n_mutually_exclusive_group,
                                 private$n_groups)
                    },
                    add_argument = function(...) {
                      private$python_code$append(sprintf("%s.add_argument(%s)", private$name,
                                                         convert_..._to_arguments("add_argument", ...)))
                      invisible(NULL)
                    },
                    add_argument_group = function(...) {
                      group_name <- paste0(private$name, "_group", private$n_groups)
                      private$n_groups <- private$n_groups + 1
                      private$python_code$append(sprintf("%s = %s.add_argument_group(%s)",
                                                         group_name, private$name,
                                                         convert_..._to_arguments("add_argument", ...)))
                      Group$new(private$python_code, group_name)
                    },
                    add_mutually_exclusive_group = function(required = FALSE) {
                      group_name <- paste0(private$name, "_mutually_exclusive_group",
                                           private$n_mutually_exclusive_groups)
                      private$n_mutually_exclusive_groups <- private$n_mutually_exclusive_groups + 1
                      private$python_code$append(sprintf("%s = %s.add_mutually_exclusive_group(%s)",
                                                         group_name, private$name,
                                                         ifelse(required, "required=True", "")))
                      Group$new(private$python_code, group_name)
                    },
                    add_subparsers = function(...) {
                      subparsers_name <- paste0(private$name, "_subparsers")
                      private$python_code$append(sprintf("%s = %s.add_subparsers(%s)",
                                                         subparsers_name, private$name,
                                                         convert_..._to_arguments("add_argument", ...)))
                      Subparsers$new(private$python_code, subparsers_name)
                    },
                    initialize = function(python_code, name,n_mutually_exclusive_groups = 0, n_groups = 0) {
                      private$python_code <- python_code
                      private$name <- name
                      private$n_mutually_exclusive_groups <- n_mutually_exclusive_groups
                      private$n_groups <- n_groups
                    }
                  ),
                  private = list(python_code = NULL, name = NULL,
                                 n_mutually_exclusive_groups = 0, n_groups = 0,n_parsers = 1)
)

parse_args_output <- function(output) {
  if (grepl("^usage:", output[1])) {
    has_positional_arguments <- any(grepl("^positional arguments:", output))
    has_optional_arguments <- any(grepl("^optional arguments:|^options:", output))
    if (has_positional_arguments || has_optional_arguments) {
      print_message_and_exit(output, "help requested:")
    } else {
      pa_stop(output, "parse error:")
    }
  } else if (grepl("^Traceback", output[1])) {
    pa_stop(output, "python error:")
  } else if (any(grepl("^SyntaxError: Non-ASCII character", output))) {
    message <- paste("Non-ASCII character detected.",
                     "If you wish to use Unicode arguments/options",
                     "please upgrade to Python 3.2+",
                     "Please see file INSTALL for more details.")
    pa_stop(message, "non-ascii character error:")
  } else if (any(grepl("^SyntaxError: positional argument follows keyword argument", output)) ||
             grepl("^SyntaxError: non-keyword arg after keyword arg", output[2])) {
    message <- "Positional argument following keyword argument."
    pa_stop(message, "syntax error:")
  } else if (grepl("^\\{|^\\[", output[1])) {
    args <- jsonlite::fromJSON(paste(output, collapse = ""))
    return(args)
  } else { # presumably version number request
    print_message_and_exit(output, "version requested:")
  }
}

# @param argument argument to be converted from R to Python
convert_argument <- function(argument, as_list = FALSE) {
  if (is.character(argument)) argument <- shQuote(argument, type = "sh")
  if (is.numeric(argument)) argument <- as.character(argument)
  if (is.logical(argument)) argument <- ifelse(argument, "True", "False")
  if (is.null(argument)) argument <- "None"
  if (as_list) {
    argument <- sprintf("[%s]", paste(argument, collapse = ", "))
  } else if (length(argument) > 1) {
    argument <- sprintf("(%s)", paste(argument, collapse = ", "))
  }
  argument
}

get_python_type <- function(type, proposed_arguments) {
  python_type <- switch(type,
                        character = "str",
                        double = "float",
                        integer = "int",
                        logical = "logical",
                        stop(paste(sprintf("type %s not supported,", type),
                                   "supported types:",
                                   "'logical', 'integer', 'double' or 'character'")))
  sprintf("type=%s", python_type)
}

should_as_list <- function(name, argument_list) {
  if (name == "default" &&
      (argument_list[["action"]] %||% "store") == "append") {
    TRUE
  } else {
    FALSE
  }
}

`%||%` <- function(x, y) if (is.null(x)) y else x # nolint

# @param mode Either "add_argument" or "ArgumentParser"
convert_..._to_arguments <- function(mode, ...) { # nolint
  
  argument_list <- list(...)
  argument_names <- names(argument_list)
  if (is.null(argument_names))
    argument_names <- rep_len("", length(argument_list))
  equals <- ifelse(argument_names == "", "", "=")
  proposed_arguments <- c()
  for (ii in seq_along(argument_list)) {
    name <- argument_names[ii]
    equal <- equals[ii]
    as_list <- should_as_list(name, argument_list)
    argument <- convert_argument(argument_list[[ii]], as_list)
    proposed_arguments <- append(proposed_arguments,
                                 paste0(name, equal, argument))
  }
  # Make sure types are what Python wants
  if (mode == "add_argument" && any(grepl("type=", proposed_arguments))) {
    ii <- grep("type=", proposed_arguments)
    type <- argument_list[[ii]]
    proposed_arguments[ii] <- get_python_type(type, proposed_arguments)
  }
  # Don't put quotes around formatter_class argument
  if (mode == "ArgumentParser" && any(grepl("formatter_class=", proposed_arguments))) {
    ii <- grep("formatter_class=", proposed_arguments)
    formatter_class <- argument_list[[ii]]
    proposed_arguments[ii] <- sprintf("formatter_class=%s", formatter_class)
  }
  # Set right default prog name if not specified, if possible
  # Do last to not screw up other fixes with prog insertion
  if (mode == "ArgumentParser" && needs_prog(argument_names)) {
    prog <- get_Rscript_filename()
    if (is.na(prog)) prog <- "PROGRAM"
    proposed_arguments <- c(sprintf("prog='%s'", prog), proposed_arguments)
  }
  return(paste(proposed_arguments, collapse = ", "))
}

needs_prog <- function(argument_names) {
  if (length(argument_names) == 0L) {
    TRUE
  } else if (argument_names[1L] == "") {
    FALSE
  } else if (any(argument_names == "prog")) {
    FALSE
  } else {
    TRUE
  }
}

# Manually copied over from getopt to eliminate it as a dependency
get_Rscript_filename <- function() { # nolint
  prog <- sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)[1])
  if (.Platform$OS.type == "windows") {
    prog <- gsub("\\\\", "\\\\\\\\", prog)
  }
  prog
}

# Internal function to check python cmd is okay
# @param python_cmd Python cmd to use
assert_python_cmd <- function(python_cmd) {
  if (!findpython::is_python_sufficient(python_cmd, required_modules = c("argparse", "json | simplejson"))) {
    stop(paste(sprintf("python executable %s either is not installed,", python_cmd),
               "is not on the path, or does not have argparse, json modules",
               "please see INSTALL file"))
  }
}

detects_python <- function() {
  python_cmd <- try(find_python_cmd())
  !inherits(python_cmd, "try-error")
}

# Internal function to find python cmd
# @param python_cmd  Python cmd to use
find_python_cmd <- function(python_cmd = NULL) {
  if (is.null(python_cmd)) {
    python_cmd <- getOption("python_cmd")
  }
  if (is.null(python_cmd)) {
    required_modules <- c("argparse", "json | simplejson")
    did_find_python3 <- findpython::can_find_python_cmd(minimum_version = "3.0",
                                                        required_modules = required_modules,
                                                        silent = TRUE)
    if (did_find_python3) {
      python_cmd <- attr(did_find_python3, "python_cmd")
    } else {
      python_cmd <- findpython::find_python_cmd(required_modules = required_modules)
    }
  }
  python_cmd
}

pa_stop <- function(message, r_note) {
  msg <- paste(c(r_note, message), collapse = "\n")
  cnd <- errorCondition(msg,
                        call = "argparse::parse_args_output(output)",
                        class = "argparse_parse_error")
  if (interactive()) {
    stop(cnd)
  } else {
    signalCondition(cnd)
    cat(message, sep = "\n", file = stderr())
    opt <- options(show.error.messages = FALSE)
    on.exit(options(opt))
    stop(cnd)
  }
}

# Internal function to print message
print_message_and_exit <- function(message, r_note) {
  if (interactive()) {
    pa_stop(message, r_note)
  } else {
    cat(message, sep = "\n")
    quit(status = 0)
  }
}

