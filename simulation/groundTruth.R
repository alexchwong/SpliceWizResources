library(fst)
library(data.table)
library(dplyr)
library(genefilter)
library(ggplot2)
library(ROCit)

library(SpliceWiz)

getGroundTruth <- function(flux_dir, ref_path) {
    pros = findSamples(flux_dir, ".pro")
    dat = fread(pros$path[1])
    expr = data.frame(transcript_id = dat$V2)
    for(i in seq_len(nrow(pros))) {
        dat = fread(pros$path[i])
        expr$TPM = dat$V5
        colnames(expr)[1 + i] = pros$sample[i]
    }
    samplenames <- colnames(expr)[-1]

    transcripts = fst::read.fst(file.path(ref_path, "fst", "Transcripts.fst"))

    missing_tid = transcripts$transcript_id[
        !(transcripts$transcript_id %in% expr$transcript_id)
    ]
    missing.df = data.frame(transcript_id = missing_tid)
    zero.mat = matrix(0, nrow = nrow(missing.df), ncol = ncol(expr) - 1)
    missing.df = cbind(missing.df, as.data.frame(zero.mat))
    colnames(missing.df)[-1] = samplenames
    expr = rbind(expr, missing.df)
    expr = expr %>% arrange(transcript_id)

    ground_truth_expr = expr
    rm(expr, missing.df, pros, transcripts, zero.mat, dat, missing_tid, i)

    ## Perform T-test on ground truth expressions, based on updated annotations
    gt_cts <- ground_truth_expr[,-1]
    rownames(gt_cts) <- ground_truth_expr$transcript_id
    gt_ttest <- test_diff_PSI(gt_cts, ref_path)

    ## Compile a list of ground truth PSIs for the rowEvents
    splice.options = fst::read.fst(file.path(ref_path, "fst", "Splice.options.fst"))
    splice.PSI = as.data.table(splice.options[, 
        c("EventType", "EventID", "isoform", "transcript_id")])
    expr = as.data.table(ground_truth_expr)
    PSI.splice = expr[splice.PSI, on = "transcript_id"]

    PSI.splice = PSI.splice[, lapply(.SD, sum, na.rm = TRUE), 
        by = c("EventType", "EventID",  "isoform"), 
        .SDcols = samplenames
    ]
    setorderv(PSI.splice, "EventID")
    mat.a = as.matrix(PSI.splice[isoform == "A", samplenames, with = FALSE])
    mat.b = as.matrix(PSI.splice[isoform == "B", samplenames, with = FALSE])
    PSI.gt = mat.a / (mat.a + mat.b)

    # Match EventName with EventID
    splice = fst::read.fst(file.path(ref_path, "fst", "Splice.fst"))
    rownames(PSI.gt) = splice$EventName[
        match(unique(PSI.splice$EventID), splice$EventID)]

    ## Calculate arithmetic mean PSI
    PSI.mean <- data.frame(
        mean_A = rowMeans(PSI.gt[,1:3]),
        mean_B = rowMeans(PSI.gt[,4:6]),
        row.names = rownames(PSI.gt)
    )
    gt_ttest$EventName <- splice$EventName[match(gt_ttest$EventID, splice$EventID)]
    gt_ttest$meanPSI_A <- PSI.mean$mean_A[match(gt_ttest$EventName, rownames(PSI.mean))]
    gt_ttest$meanPSI_B <- PSI.mean$mean_B[match(gt_ttest$EventName, rownames(PSI.mean))]
    gt_ttest$deltaPSI <- gt_ttest$meanPSI_A - gt_ttest$meanPSI_B
    
    gt <- list(
        gt_expr = ground_truth_expr,
        gt_PSI = PSI.gt,
        gt_diff = gt_ttest
    )
    return(gt)
}

generatePSIerror <- function(
        se, DE_results, gt,
        samplenames = colnames(se),
        reference_path
) {
    res <- DE_results %>% filter(EventType != "IR")
    splice <- fst::read.fst(file.path(reference_path, "fst", "Splice.fst"))
    PSI <- makeMatrix(se, res$EventName, 
        logit_max = 10, na.percent.max = 1.0)
    # PSI <- PSI[rownames(PSI) %in% rownames(gt_PSI),]
    PSI_mean <- rowMeans(PSI[,samplenames])
    PSI.gt.compare <- rowMeans(gt$gt_PSI[rownames(PSI), samplenames])
    PSI.diff <- PSI.gt.compare - PSI_mean
    PSIerror <- data.frame(
        EventName = rownames(PSI),
        mean_diff = abs(PSI.diff),
        splice_type = splice$EventType[match(rownames(PSI), splice$EventName)]
    ) %>% arrange(mean_diff) %>% filter(is.finite(mean_diff))
    return(PSIerror)
}

plotPSIerror <- function(PSIerror) {
    SW.meandiff.summa <- PSIerror %>% group_by(splice_type) %>%
        summarize(
            mean_error = mean_diff, 
            cumfreq = seq_len(n()) * 100 / n()
        ) %>% ungroup() %>% arrange(mean_error)
    p <- ggplot(SW.meandiff.summa, 
        aes(x = mean_error * 100, y = cumfreq, color = splice_type)) + 
        geom_line() + theme_white_legend
    return(p)
}

getPSIerrorAUC <- function(PSIerror, label = "Unfiltered") {
    md <- rbind(PSIerror, PSIerror %>% mutate(splice_type = "All Events"))
    summa <- md %>% 
        dplyr::select(splice_type, mean_diff) %>% 
        summa_ROC("splice_type") %>%
        mutate(method = label) %>% 
        clean_ROC_summa(splice_type, method)
    PSIerrorAUC <- get_ROC_values(summa)
    return(PSIerrorAUC)
}

generateScoresAndClass <- function(
        DE_results, gt, reference_path,
        padj_threshold = 0.05,
        log2_min_threshold = -1,
        deltaPSI_threshold = 0.05
) {
    splice <- fst::read.fst(file.path(reference_path, "fst", "Splice.fst"))
    gt_ttest <- gt$gt_diff
    gt_ttest$log2_min <- with(gt_ttest, 
        ifelse(log2_mean_A < log2_mean_B, log2_mean_A, log2_mean_B)
    )
    gt_ttest$evaluate <- with(gt_ttest, log2_min > log2_min_threshold)
    gt_ttest$truth <- with(gt_ttest, 
        padj < padj_threshold & 
        log2_min > log2_min_threshold & 
        abs(deltaPSI) > deltaPSI_threshold
    )
    gt_ttest$EventName <- splice$EventName[match(gt_ttest$EventID, splice$EventID)]
    
    DE_results <- DE_results %>% filter(EventType != "IR")
    DE_results$truth <- gt_ttest$truth[match(DE_results$EventName, gt_ttest$EventName)]
    DE_results$evaluate <- gt_ttest$evaluate[match(DE_results$EventName, gt_ttest$EventName)]
    DE_results$EventID <- gt_ttest$EventID[match(DE_results$EventName, gt_ttest$EventName)]
    
    if("pvalue" %in% colnames(DE_results))
        DE_results <- DE_results %>% dplyr::rename(P.Value = pvalue)
    
    ScoresAndClass <- DE_results[, c("EventID", "P.Value", "truth", "evaluate")] %>%
        mutate(splice_type = tstrsplit(EventID, split="#", fixed=TRUE)[[1]])
    ScoresAndClass$P.Value[!is.finite(ScoresAndClass$P.Value)] <- 1
    ScoresAndClass <- ScoresAndClass %>% arrange(P.Value) %>% 
        filter(evaluate) %>% dplyr::select(splice_type, truth, P.Value)
    return(ScoresAndClass)
}

generateROCdata <- function(ScoresAndClass) {
    # Overall
    roc <- suppressWarnings(
        ROCit::rocit(-ScoresAndClass$P.Value, ScoresAndClass$truth)
    )
    roc_df <- rocit_2_df(roc) %>% mutate(splice_type = "All Events")
    ROCdata <- roc_df
    # Modality Specific
    for(st in unique(ScoresAndClass$splice_type)) {
        sc_subset <- ScoresAndClass %>% filter(splice_type == st)
        roc <- suppressWarnings(
            ROCit::rocit(-sc_subset$P.Value, sc_subset$truth)
        )
        roc_df <- rocit_2_df(roc) %>% mutate(splice_type = st)
        ROCdata <- rbind(ROCdata, roc_df)
    }
    return(ROCdata)
}

plotROCdata <- function(ROCdata) {
    p <- ggplot(ROCdata, aes(y = TPR, x = FPR)) + 
        geom_line() + theme_white_legend +
        facet_wrap(vars(splice_type))
    p
}

getAUROC <- function(ROCdata) {
    ret <- ROCdata %>% group_by(splice_type) %>% 
        summarize(AUROC = max(cum_AUC)) %>% ungroup()
    return(as.data.frame(ret))
}

## Internal functions

# Perform a t-test for a given set of expressions
test_diff_PSI <- function(cts, reference_path) {

    splice.options = fst::read.fst(file.path(reference_path, "fst/Splice.options.fst"))
    junctions = fst::read.fst(file.path(reference_path, "fst/junctions.fst"))
    exons = fst::read.fst(file.path(reference_path, "fst/Exons.fst"))

    n_reps = ncol(cts) / 2
    # Convert these back to TPMs
    # cts.TPM = cts / rowMeans(cts) * mean_TPM
    cts.TPM = t(t(cts) / colSums(cts)) * 1e6
    cts.TPM[!is.finite(cts.TPM)] = 0

    # Simulate ASE using NxtIRF table
    SO.expr = as.data.table(
        splice.options[, c("EventType", "EventID", "transcript_id", "isoform")])

    cts.TPM.DT = as.data.table(cts.TPM)
    cts.TPM.DT$transcript_id = rownames(cts.TPM)

    temp = cts.TPM.DT[SO.expr, on = "transcript_id"]

    temp.mat = as.matrix(temp[, seq_len(n_reps * 2), with = FALSE])
    temp.mat[is.na(temp.mat)] = 0
    temp.mat = as.data.table(temp.mat)

    SO.expr = cbind(SO.expr, temp.mat)

    SO.expr.summa = copy(SO.expr)
    SO.expr.summa$transcript_id = NULL
    SO.expr.summa = SO.expr.summa[, lapply(.SD, sum, na.rm = TRUE), 
        by = c("EventType", "EventID",  "isoform"), 
        .SDcols = colnames(cts)
    ]
    SO.expr.summa = unique(SO.expr.summa, by = c("EventID", "isoform"))
    setorder(SO.expr.summa, "EventID")
    SO.expr.A = SO.expr.summa[isoform == "A"]
    SO.expr.B = SO.expr.summa[isoform == "B"]

    mat.a = as.matrix(SO.expr.A[, colnames(cts), with = FALSE])
    mat.b = as.matrix(SO.expr.B[, colnames(cts), with = FALSE])
    psi = mat.a / (mat.a + mat.b)

    ### Survey how many ASEs are already differentially expressed:

    # Cap at logit_psi < 10
    psi_logit = qlogis(psi)
    psi_logit[psi_logit > 10] = 10
    psi_logit[psi_logit < -10] = -10

    # T-test of logit_psi
    t_test = genefilter::rowttests(
        psi_logit, factor(rep(c("A", "B"), each = n_reps)), na.rm = TRUE)
    t_test = as.data.frame(t_test)

    t_test$log2_mean_A = log2(rowMeans(mat.a) + 0.01)
    t_test$log2_mean_B = log2(rowMeans(mat.b) + 0.01)
    t_test$log2_mean = log2(rowMeans(mat.b) + rowMeans(mat.a) + 0.01)
    t_test$mean_PSI = with(t_test, rowMeans(mat.a) / (rowMeans(mat.a) + rowMeans(mat.b)))
    sum(t_test$p.value < 0.05, na.rm=TRUE)

    SO.expr.final = SO.expr.A[,1:2]
    SO.expr.final$mean_TPM = rowMeans(mat.b) + rowMeans(mat.a)
    SO.expr.final$mean_TPM_A = rowMeans(mat.a)
    SO.expr.final$mean_TPM_B = rowMeans(mat.b)

    t_test = cbind(t_test, SO.expr.final )
    t_test.ranked = t_test %>% arrange(p.value)
    t_test.ranked$padj = p.adjust(t_test.ranked$p.value, method = "BH")
    
    return(t_test.ranked)
}

summa_ROC <- function(df, bycol = "splice_type") {
    df$group_col <- "same"

    if(bycol != "") df$group_col <- df[, bycol]
    if(!("method" %in% names(df))) df$method <- "none"

    df <- df %>% filter(is.finite(mean_diff)) %>%
        arrange(mean_diff) %>% 
        group_by(group_col, method) %>%
        summarize(
            mean_error = mean_diff, cumfreq = seq_len(n()) * 100 / n()
        ) %>% ungroup()
    colnames(df)[colnames(df) == "group_col"] = bycol
    return(df)
}

clean_ROC_summa <- function(df, ...) {
    df.spawn <- df %>%
        dplyr::select( ...) %>% 
        distinct(..., .keep_all = TRUE) %>%
        mutate(mean_error = 1, cumfreq = 100)
    df <- rbind(df, df.spawn)
    return(df)
}

get_ROC_values <- function(summa) {
  # method, splice_type
    mets <- unique(summa$method)
    classes <- unique(summa$splice_type)
    mat <- matrix(nrow = length(classes), ncol = length(mets))

    for(i in seq_len(length(classes))) {
        for(j in seq_len(length(mets))) {
            met <- mets[j]
            as <- classes[i]
            if(any(met %in% (summa %>% filter(splice_type == as))$method)) {
                mat[i, j] <- get_AUC(
                    summa %>% filter(splice_type == as & method == met)
                )
            } 
        }
    }
    rownames(mat) <- classes
    colnames(mat) <- mets
    mat
}

get_AUC <- function(summa) {
    # Assume summa is filtered by method and splice_type
    summa$mean_error[summa$mean_error > 100] <- 100
    summa <- summa %>% arrange(mean_error) %>%
        mutate(last_mean_error = lag(mean_error, default = 0)) %>%
        mutate(last_cumfreq = lag(cumfreq, default = 0))
    summa$AUC <- with(summa, 
        (cumfreq + last_cumfreq) * (mean_error - last_mean_error) / 2)
    summa$cum_AUC <- cumsum(summa$AUC)
    return(max(summa$cum_AUC) / 100)
}

rocit_2_df <- function(rocit) {
    df <- data.frame(
        TPR = rocit$TPR,
        FPR = rocit$FPR,
        cutoff = -rocit$Cutoff
    )
    df <- df %>% arrange(FPR) %>%
        mutate(last_FPR = lag(FPR, default = 0)) %>%
        mutate(last_TPR = lag(TPR, default = 0)) %>%
        mutate(AUC = (last_TPR + TPR) * (FPR - last_FPR) / 2)
    df$cum_AUC <- cumsum(df$AUC)
    df
}