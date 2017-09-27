#' Calculates leaf traits for a directory of leaf scans
#' either made with a white or black background.
#'
#' @param path: path of image files to process (should only contain images),
#' default is your home directory
#' @param dpi: dots per inch resolution of the scan (default = 300)
#' @param background: background matte used ("white" or "black")
#' @param reference: is there a reference strip to the left which
#' should be excluded, provide a number of pixels to do so (default = 0)
#' @param out_path: output directory where to save the data if not returned
#' to an R variable
#' @param plot: plot the processed images for reference (default = FALSE)
#' @return A nested list (data structure) listing the extracting traits together
#' with the filename of the processed file.
#' @keywords leaf traits, specific leaf area, leaf area, Gcc, Rcc
#' @export
#' @examples
#'
#' \dontrun{
#'
#' leaf_statistics = calculate_leaf_traits(path = "/mydata")
#'
#' }

calculate_leaf_traits = function(path = "~/Desktop/testdir/",
                                 dpi = 300,
                                 background = "white",
                                 reference = 0,
                                 out_path = NULL,
                                 plot = TRUE){

  # routine to calculate Gcc / Rcc
  getCover = function(img, roi, updatevalue){

    # read in the image and ROI mask
    m = roi == updatevalue
    m[m == 0] = NA

    # read image bands
    r = raster(img, band=1)
    g = raster(img, band=2)
    b = raster(img, band=3)

    # mask all values which are of no
    # importance
    r_s = mask(r, m, maskevalue = NA)
    g_s = mask(g, m, updatevalue = NA)
    b_s = mask(b, m, updatevalue = NA)

    # average across the ROI
    r_m = cellStats(r_s, stat = mean)
    g_m = cellStats(g_s, stat = mean)
    b_m = cellStats(b_s, stat = mean)

    # calculate the gcc (% green)
    gcc =  g_m/(r_m + g_m + b_m)
    rcc =  r_m/(r_m + g_m + b_m)

    # combine raw and calculated numbers
    data = as.numeric(list(r_m,g_m,b_m,gcc,rcc))
    return(data)
  }

  # list all files in the directory
  # quit if no files are detected
  files = list.files(path, "*", full.names = TRUE)

  if(length(files) == 0){
    stop("No files in specified directory!")
  }

  batch_output = lapply(files,function(file){

    # read in image
    im = load.image(file)

    if (reference != 0 || is.null(reference)){
      px = Xc(im) <= reference
      # blank out reference panels
      # could be done manually as well
      if (background == "white"){
        bg = imfill(dim=dim(im),val=c(1,1,1))
        msk = as.cimg(px)
        im = bg*msk+(1-msk)*im
      }else{
        im[px] = 0
      }
    }

    # convert to CIELAB colour space,
    # then create a data.frame with three colour channels as columns
    d = sRGBtoLab(im) %>% as.data.frame(wide="c")%>%
      dplyr::select(-x,-y)
    L = channel(RGBtoHSL(im),3)
    # convert to raster format
    L = raster(t(as.matrix(L)))
    extent(L) = extent(0,ncol(L),0,nrow(L))

    # run k-means with 2 centers
    km = kmeans(d,2)

    # turn cluster index into an image
    seg = as.cimg(km$cluster,dim=c(dim(im)[1:2],1,1))

    # convert to raster format
    seg = raster(t(as.matrix(seg)))
    extent(seg) = extent(0,ncol(seg),0,nrow(seg))

    # calculate zonal statistics
    zstats = as.data.frame(zonal(L, seg, fun = "mean", na.rm = TRUE))

    if (background == "white"){
      seg[seg == which(zstats$mean == max(zstats$mean, na.rm = TRUE))] = NA
    } else {
      seg[seg == which(zstats$mean == min(zstats$mean, na.rm = TRUE))] = NA
    }

    # clump the data
    # kick out anything smaller than 100 px
    c = clump(seg, gaps = FALSE)
    extent(c) = extent(0,ncol(c),0,nrow(c))
    for(i in which(freq(c)[,2] < 100) ){
     c[c == i] = NA
    }

    # get centroids
    clump_id = getValues(c)
    xy = xyFromCell(c,1:ncell(c))
    df = na.omit(data.frame(xy, clump_id))
    centroids = by(df,INDICES = df$clump_id, function(xy){
      x = median(xy$x[1])
      y = median(xy$y[2])
      data.frame(x,y)}
       )
    centroids = do.call("rbind",centroids)

    # discard noise
    pixel_counts = as.data.frame(freq(c, useNA = 'no'))
    results = lapply(pixel_counts$value,function(updatevalue){
                getCover(file, c, updatevalue)
              })

    # bind things
    pixel_stats = data.frame(basename(file),
                        centroids,
                        pixel_counts$count,
                        do.call("rbind", results))
    colnames(pixel_stats) = c("filename","x","y","count","r","g","b","gcc","rcc")

    # convert counts to mm^2
    pixel_stats$square_cm = pixel_stats$count / (dpi / 2.54)^2

    # plot results for feedback
    if (plot){
      plotRGB(brick(file))
      text(pixel_stats$x,
           pixel_stats$y,
           adj = c(0,3),
           col = "red")

      text(pixel_stats$x,
           pixel_stats$y,
           paste("Gcc: ", round(pixel_stats$gcc,3)),
           adj = c(0,5),
           col = "red")

      text(pixel_stats$x,
           pixel_stats$y,
           paste("mm^2: ", round(pixel_stats$square_cm,2)),
           adj = c(0,7),
           col = "red")
    }
    return(pixel_stats)
  })

  # output results to file or to R workspace
  if (is.null(out_path)){
    return(batch_output)
  } else {
    saveRDS(batch_output, sprintf("%s/leaf_trait_image_analysis.rds", out_path))
  }
}
