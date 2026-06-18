#' @title Generate a Rule-Based Breeding Matrix
#' @description This function simulates a breeding matrix based on a set of
#'   rules and a Poisson distribution for the number of mates per individual.
#'   This method is known to introduce a systematic bias due to the order of
#'   operations and is used here to demonstrate the superiority of a more
#'   rigorous approach.
#' @param moms The number of female parents in the population.
#' @param dads The number of male parents in the population.
#' @param lambda.low The lower bound of the lambda value for the Poisson
#'   distribution used to determine the number of mates.
#' @param lambda.high The upper bound of the lambda value for the Poisson
#'   distribution used to determine the number of mates.
#' @return A binary breeding matrix with `dads` rows and `moms` columns, where
#'   a '1' indicates a mating occurred.
#' @examples
#' # Generate a simple breeding matrix
#' my_matrix <- brd.mat(moms = 50, dads = 50, lambda.low = 1, lambda.high = 1)
#' # Note: This is for demonstration of the old method's internal biases.
brd.mat <- function(moms = 10, dads = 10, lambda.low = 1, lambda.high = 1) {
  mat <- matrix(data = 0, nrow = dads, ncol = moms)
  lambda.mates <- sample(lambda.low:lambda.high, size = 1)
  dad.mate.probs <- data.frame(mates = 1:dads)
  suppressWarnings(dad.mate.probs$prob <- (lambda.mates^dad.mate.probs$mates) * (exp(-lambda.mates)) / factorial(dad.mate.probs$mates))
  dad.mate.probs <- dad.mate.probs[dad.mate.probs$prob != 0 & !is.na(dad.mate.probs$prob),]
  mom.mate.probs <- data.frame(mates = 1:moms)
  suppressWarnings(mom.mate.probs$prob <- (lambda.mates^mom.mate.probs$mates) * (exp(-lambda.mates)) / factorial(mom.mate.probs$mates))
  mom.mate.probs <- mom.mate.probs[mom.mate.probs$prob != 0 & !is.na(mom.mate.probs$prob),]
  mom.order <- paste(1:moms, "Female", sep = "_")
  dad.order <- paste(1:dads, "Male", sep = "_")
  parent.order <- sample(x = c(mom.order, dad.order), size = moms + dads, replace = FALSE)
  tmp <- data.frame(parent = parent.order, stringsAsFactors = F)
  tmp$sex <- ifelse(grepl(pattern = "Female", x = tmp$parent), "Female", "Male")
  tmp$number <- as.numeric(gsub(pattern = "_.*", replacement = "", x = tmp$parent))
  tmp$mates_before <- NA
  for (i in 1:nrow(tmp)) {
    if (grepl(pattern = "Female", x = tmp$parent[i])) {
      my.col <- tmp$number[i]
      mates <- sample(dad.mate.probs$mates, size = 1, prob = dad.mate.probs$prob)
      tmp$mates_before[i] <- mates
      current.mates.list <- which(mat[, my.col] == 1)
      if (length(current.mates.list) == 0) {
        if (mates > 0) {
          my.rows <- sample(x = dad.mate.probs$mates, size = mates, replace = F)
          mat[my.rows, my.col] <- 1
        }
      } else {
        if (length(current.mates.list) < mates) {
          pdm2 <- dad.mate.probs$mates[!(dad.mate.probs$mates %in% current.mates.list)]
          mates2 <- mates - length(current.mates.list)
          if(length(pdm2) > 0 && mates2 > 0) {
            my.rows <- sample(x = pdm2, size = min(mates2, length(pdm2)))
            mat[my.rows, my.col] <- 1
          }
        }
        if (length(current.mates.list) > mates) {
          cml2 <- sample(x = current.mates.list, size = mates)
          rem.mates <- current.mates.list[!(current.mates.list %in% cml2)]
          mat[rem.mates, my.col] <- 0
        }
      }
    } else {
      my.row <- tmp$number[i]
      mates <- sample(mom.mate.probs$mates, size = 1, prob = mom.mate.probs$prob)
      tmp$mates_before[i] <- mates
      current.mates.list <- which(mat[my.row, ] == 1)
      if (length(current.mates.list) == 0) {
        if (mates > 0) {
          my.cols <- sample(x = mom.mate.probs$mates, size = mates, replace = F)
          mat[my.row, my.cols] <- 1
        }
      } else {
        if (length(current.mates.list) < mates) {
          pmm2 <- mom.mate.probs$mates[!(mom.mate.probs$mates %in% current.mates.list)]
          mates2 <- mates - length(current.mates.list)
          if(length(pmm2) > 0 && mates2 > 0) {
            my.cols <- sample(x = pmm2, size = min(mates2, length(pmm2)))
            mat[my.row, my.cols] <- 1
          }
        }
        if (length(current.mates.list) > mates) {
          cml2 <- sample(x = current.mates.list, size = mates)
          rem.mates <- current.mates.list[!(current.mates.list %in% cml2)]
          mat[my.row, rem.mates] <- 0
        }
      }
    }
  }
  tmp <- tmp[order(tmp$number),]
  tmp <- tmp[order(tmp$sex),]
  tmp$order <- 1:nrow(tmp)
  tmp$mates_after <- c(colSums(mat), rowSums(mat))
  tmp$diff <- tmp$mates_after - tmp$mates_before
  tmp <- tmp[order(tmp$diff, decreasing = T),]
  i <- 0
  black.list <- NULL
  blacklist.carryover <- NULL
  while(sum(tmp$diff)>0 && i < dads + moms) {
    i <- i + 1
    tmp1 <- tmp[!(tmp$parent %in% black.list), ]
    if(length(blacklist.carryover) > 0) {
      tmp1 <- tmp1[!(tmp1$parent %in% blacklist.carryover), ]
    }
    if(nrow(tmp1) == 0 || tmp1[1,]$diff < 0) {
      break
    }
    tmp1 <- tmp1[1, ]
    if(tmp1$sex == "Female") {
      if(tmp1$diff > 0) {
        my.males <- which(mat[, tmp1$number] > 0)
        inflated <- tmp[tmp$sex == "Male" & (tmp$number %in% my.males) & tmp$diff > 0, ]
        if(nrow(inflated) == tmp1$diff) {
          mat[inflated$number, tmp1$number] <- 0
        } else if(nrow(inflated) > tmp1$diff) {
          my.removals <- sample(x = inflated$number, size = tmp1$diff)
          mat[my.removals, tmp1$number] <- 0
        } else if(nrow(inflated) < tmp1$diff && nrow(inflated) > 0) {
          mat[inflated$number, tmp1$number] <- 0
        } else if(nrow(inflated) == 0) {
          blacklist.carryover <- c(blacklist.carryover, tmp1$parent)
        }
      }
    } else {
      if(tmp1$diff > 0) {
        my.females <- which(mat[tmp1$number,] > 0)
        inflated <- tmp[tmp$sex == "Female" & (tmp$number %in% my.females) & tmp$diff > 0, ]
        if(nrow(inflated) == tmp1$diff) {
          mat[tmp1$number, inflated$number] <- 0
        } else if(nrow(inflated) > tmp1$diff) {
          my.removals <- sample(x = inflated$number, size = tmp1$diff)
          mat[tmp1$number, my.removals] <- 0
        } else if(nrow(inflated) < tmp1$diff && nrow(inflated) > 0) {
          mat[tmp1$number, inflated$number] <- 0
        } else if(nrow(inflated) == 0) {
          blacklist.carryover <- c(blacklist.carryover, tmp1$parent)
        }
      }
    }
    tmp$mates_after <- c(colSums(mat), rowSums(mat))
    tmp$diff <- tmp$mates_after - tmp$mates_before
    tmp <- tmp[order(tmp$diff, decreasing = T), ]
    black.list <- tmp$parent[tmp$diff == 0]
  }
  return(mat)
}
