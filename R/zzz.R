# .setup.pp.data <- function(data, responsename, pp) {
# 
# nodes <- pp$nodes
# ny <- pp$ny
# weights <- pp$weights
# threshold <- pp$threshold
# r <- pp$r
# 
# if (is.null(nodes) & is.null(ny))
#   stop("Cannot have NULL pp.args$ny and pp.args$nodes.")
# if (is.null(threshold) & is.null(r)) 
#   stop("Both pp$threshold and pp$r cannot be NULL")
# 
# if (!is.null(pp$nodes)) { # start stuff with nodes
# 
# ## option for integration over nodes added to evgam_0.1.2
# # otherwise integration over pp$id, and assumed constant
# print("Testing pp$nodes")
# 
# if (!inherits(nodes, "data.frame")) {
# 
# print("Nodes are a list")
# 
# # Assume nodes for each variable given, and then use expand.grid
# # to get all combinations
# # Second column of each list element taken as weights, if present
# 
# if (!inherits(nodes, "list")) {
#   stop('Unrecognized class(pp$nodes)[1] not in c("data.frame", "list")')
# } else {
#   nodes <- lapply(nodes, as.matrix)
#   for (i in seq_along(nodes)) if (ncol(nodes[[i]]) == 1) nodes[[i]] <- cbind(nodes[[i]], 1)
#   wts <- lapply(nodes, function(x) x[,2])
#   wts <- Reduce("*", expand.grid(wts))
#   nodes <- lapply(nodes, function(x) x[,1])
#   nodes <- expand.grid(nodes)
#   nd <- nrow(nodes)
#   if (is.null(weights)) {
#     weights <- 1
#   } else {
#     if (length(weights) > 1) {
#       stop('Cannot have class(pp$nodes)[1] == "list" and length(pp$weights) > 1')
#     } else {
#       wts <- weights * wts
#     }
#   }
# }
# 
# # Now deal with threshold
# 
# if (!is.null(threshold)) {
#   if (length(threshold) > 1) {
#     stop('Cannot have length(threshold) > 1')
#   } else {
#     if (inherits(threshold, c("numeric", "integer"))) {
#       data <- subset(data, data[,responsename] >= threshold)
#       nodes[,responsename] <- threshold
#     } else {
#       if (inherits(threshold, "character")) {
#         if (threshold %in% names(nodes)) {
#           names(nodes)[names(nodes) == threshold] <- responsename
#         } else {
#           if (threshold %in% names(data)) {
#             data <- subset(data, data[,responsename] >= data[,threshold])
#           } else {
#             stop(paste(threshold, "not in names(data)"))
#           }
#           stop(paste(threshold, "not in names(nodes)"))
#         }
#       } else {
#         stop('class(pp$threshold) not in c("numeric", "character")')
#       }
#     }
#   }    
# } else {
#   if (!is.null(pp$r)) {
#     threshold <- sort(data[,responsename], decreasing=TRUE)[pp$r]
#     nodes[,responsename] <- threshold
#     data <- subset(data, data[,responsename] >= threshold)
#   } else {
#     stop('Both pp$r and pp$threshold cannot be NULL') # covered at start
#   }
# }
#   
# #   class(threshold)[1] == "character") {
# #     if (length(threshold) == 1) {
# #       if (threshold %in% names(nodes)) {
# #         if (threshold %in% names(data)) {
# #           data <- subset(data, data[,responsename] >= data[,threshold])
# #         } else {
# #           stop(paste(threshold, "not in names(data)"))
# #         }
# #       } else {
# #         stop(paste(threshold, "not in names(nodes)"))
# #       }
# #     } else {
# #       stop('class(threshold)[1] == "character" but length(threshold) != 1')
# #     }
# #   } else {
# #     if (class(pp$threshold)[1] %in% c("numeric", "integer")) {
# #       if (length(pp$threshold) == 1) {
# #         nodes[,responsename] <- pp$threshold   
# #       } else {
# #         stop('Cannot have class(nodes)[1] == "list" and length(pp$threshold) > 1')
# #       }
# #     } else {
# #       stop('class(pp$threshold)[1] not in c("numeric", "integer")')
# #     }
# #   }
# # } else {
# #   if (!is.null(pp$r)) {
# #     nodes[,responsename] <- sort(data[,responsename], decreasing=TRUE)[pp$r]
# #   } else {
# #     stop('Both pp$r and pp$threshold cannot be NULL')
# #   }
# # }
# 
# } else {
# 
# print("Nodes are a data.frame")
# 
# # Assumes all nodes given
# 
# nd <- nrow(nodes)
# 
# # make integration weights
# # (this might become a bit more sophisticated, like ppgam, in future)
# if (is.null(pp$weights)) {
#   wts <- rep(1, nrow(nodes))
# } else {
#   if (length(pp$weights) == 1 & inherits(pp$weights, "character")) {
#     wts <- nodes[,pp$weights]
#     nodes[,pp$weights] <- NULL
#   } else {
#     if (length(pp$weights) == 1 & inherits(pp$weights, c("numeric", "integer"))) {
#       wts <- rep(pp$weights, nd)
#     } else {
#       if (length(pp$weights) > 1 & inherits(pp$weights, c("numeric", "integer"))) {
#         if (length(pp$weights) != nd) {
#           stop("nrow(pp$nodes) and length(pp$weights) not compatible")
#         } else {
#           wts <- pp$weights
#         }
#       }
#     }
#   }
# }
# 
# # Now deal with threshold
# 
# if (!is.null(threshold)) {
#   if (length(threshold) > 1) {
#     stop('Cannot have length(threshold) > 1')
#   } else {
#     if (inherits(threshold, c("numeric", "integer"))) {
#       nodes[,responsename] <- threshold
#       data <- subset(data, data[,responsename] >= threshold)
#     } else {
#       if (inherits(threshold, "character")) {
#         if (threshold %in% names(nodes)) {
#           names(nodes)[names(nodes) == threshold] <- responsename
#         } else {
#           if (threshold %in% names(data)) {
#             data <- subset(data, data[,responsename] >= data[,threshold])
#           } else {
#             stop(paste(threshold, "not in names(data)"))
#           }
#           stop(paste(threshold, "not in names(nodes)"))
#         }
#       } else {
#         stop('class(pp$threshold) not in c("numeric", "character")')
#       }
#     }
#   }    
# } else {
#   if (!is.null(pp$r)) {
#     threshold <- sort(data[,responsename], decreasing=TRUE)[pp$r]
#     nodes[,responsename] <- threshold
#     data <- subset(data, data[,responsename] >= threshold)
#   } else {
#     stop('Both pp$r and pp$threshold cannot be NULL') # covered at start
#   }
# }
#   
# # if (!(responsename %in% names(nodes))) {
# #   if (class(pp$threshold)[1] %in% c("numeric", "integer")) {
# #     if (length(pp$threshold) == 1) {
# #       nodes[,responsename] <- pp$threshold   
# #     } else {
# #       if (length(pp$threshold) != nd) {
# #         stop('Need length(pp$threshold) in c(1, nrow(nodes)) for class(nodes)[1] == "data.frame".')
# #       } else {
# #         nodes[,responsename] <- pp$threshold   
# #       }
# #     }
# #   } else {
# #     if (!is.null(pp$r)) {
# #       nodes[,responsename] <- sort(data[,responsename], decreasing=TRUE)[pp$r]
# #       data <- subset(data, data[,responsename] >= nodes[1, responsename])
# #     } else {
# #       stop('Neither pp$r nor pp$threshold acceptably provided')
# #     }
# #   }
# # }
# 
# }
# 
# out <- data
# attr(out, "weights") <- wts
# attr(out, "quad") <- nodes
# if (is.null(pp$cens)) {
#   attr(out, "cens") <- NULL
# } else {
#   attr(out, "cens") <- data[out$row, pp$cens]
# }
# # if (is.null(pp$weights)) {
#   attr(out, "cweights") <- rep(1, nrow(out))
# # } else {
# # #   attr(out, "cweights") <- pp$weights[out$row]
# #   attr(out, "cweights") <- rep(pp$weights, nrow(out))
# # }
# if (!is.null(pp$exi)) attr(out, "exi") <- data[,pp$exi]
# print(attributes(out)[c("weights", "quad", "cens")])
# 
# # last line of stuff with nodes
# 
# } else {
# ## simple constant partial point process
# data$row <- seq_len(nrow(data))
# ds <- split(data, data[,pp$id])
# wts <- pp$ny
# if (length(wts) == 1) {
#   wts <- rep(wts, length(ds))
# } else {
#   wts <- wts[match(names(ds), names(wts))]
# }
# nobs2 <- sapply(ds, nrow)
# ## start of original r-largest order statistic stuff
# data.quad <- do.call(rbind, lapply(ds, function(x) x[1,]))
# if (!is.null(pp$r)) {
# enough <- nobs2 > pp$r
# if (any(!enough)) warning(paste(sum(!enough), "unique pp.args$id removed for having fewer than r observations."))
# ds <- ds[enough]
# wts <- wts[enough]
# nid <- sum(enough)
# data.quad <- data.quad[enough,]
# if (pp$r != -1) {
#     du <- sapply(ds, function(x) x[order(x[,responsename], decreasing=TRUE)[pp$r], responsename])
# } else {
#     du <- sapply(ds, function(x) min(x[, responsename]))
# }
# ## end of ...
# } else {
# ## start of specified threshold stuff
# du <- sapply(ds, function(x) x[1, pp$threshold])
# nid <- length(du)
# ## end of specified threshold stuff
# }
# data.quad[,responsename] <- du
# ds <- lapply(seq_len(nid), function(i) subset(ds[[i]], ds[[i]][,responsename] >= du[i]))
# out <- dfbind(ds)
# attr(out, "weights") <- wts
# attr(out, "quad") <- data.quad
# if (is.null(pp$cens)) {
#   attr(out, "cens") <- NULL
# } else {
#   attr(out, "cens") <- data[out$row, pp$cens]
# }
# if (is.null(pp$weights)) {
#   attr(out, "cweights") <-rep(1, length(out$row))
# } else {
#   attr(out, "cweights") <- pp$weights[out$row]
# }
# }
# out
# }
