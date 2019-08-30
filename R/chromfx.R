############################################################################
## Ashley Stephen Doane
## Weill Cornell Medicine  asd2007@med.cornell.edu

## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################





#' @name readPeakSummits
#' @title readPeakSummits
#' @description
#' function that reads ATAC peaks summits
#' @keywords peaks
#' @importFrom rtracklayer import.bed
#' @importFrom rtracklayer import.bed
#' @importFrom GenomeInfoDb seqlevels
#' @param psum paths to sample peak summits bed file
#' @param genome reference genome, either "hg3* (default) or "mm10"
#' @return GRanges
#' @export
readPeakSummits <- function(psum, genome="hg38"){
    if (genome=="mm10"){
        slevels <- mm10s
    } else {
        slevels <- hg38s
    }
    summG = list()
    for (p in names(psum)){
        print(p)
        z  = rtracklayer::import.bed(psum[[p]])
        GenomeInfoDb::seqlevels(z, pruning.mode="coarse") <- slevels
        zz <- GenomicRanges::resize(z, width=500, fix="center")
        zz$scoreO <- zz$score
        zz$score <- zz$score / length(zz)  ## normlaiza seq depth
        names(zz) <- basename(zz$name)
        summG[[p]] <- trim(zz) ## per sample unique summits
    }
    unlist(summG)
}


keepone <- function(gr, hitlist, FUN=which.max) {
    idx0 <- as(FUN(extractList(gr$score, hitlist)), "List")
    idx1 <- unlist(extractList(seq_along(gr), hitlist)[idx0])
    ## FIXME: what about NA's when there are no matching ranges?
    grx =  gr[idx1]
    grx[unique(names(grx))]
}





#' @name getAtlasPeaks
#' @title getHic
#' @description
#' function that takes a GRanges of peak summits and geenerates fixed width peak atlas
#' @keywords straw
#' @param peaksummits GRanges of peak summits with coloumn "score"
#' @return GRanges
#' @export
#' @author Ashley S Doane
getAtlasPeaks <- function(peaksummits)
{
    keepone <- function(gr, hitlist, FUN=which.max) {
        idx0 <- as(FUN(extractList(gr$score, hitlist)), "List")
        idx1 <- unlist(extractList(seq_along(gr), hitlist)[idx0])
        ## FIXME: what about NA's when there are no matching ranges?
        grx =  gr[idx1]
        grn <- grx[unique(names(grx))]
        return(grn)
    }

    gr = peaksummits[order(peaksummits$score, decreasing=TRUE)]
    #gr <- resize(gr, fix="center", width=500)
    grr <- GenomicRanges::reduce(gr)
    r1 = length(gr)
    gr = gr[order(gr$score, decreasing=TRUE)]
    hitlist <- as(GenomicRanges::findOverlaps(gr), "List")
    gr <- keepone(gr, hitlist)
    r2 <- length(gr) - 1
    while (r1 != r2) {
        r1 = length(gr)
        gr = gr[order(gr$score, decreasing=TRUE)]
        hitlist <- as(findOverlaps(gr), "List")
        gr <- keepone(gr, hitlist)
        r2 = length(gr)
        print(c(r1, r2))
    }
    return(gr)
}