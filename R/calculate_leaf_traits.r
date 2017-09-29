#' Calculates leaf traits for a directory of leaf scans
#' either made with a white or black background.
#'
#' @param path path of image files to process (should only contain images),
#' default is your home directory
#' @param dpi dots per inch resolution of the scan (default = 300)
#' @param background background matte used ("white" or "black")
#' @param reference is there a reference strip to the left which
#' should be excluded, provide a number of pixels to do so (default = NULL)
#' @param min_size minimum size of a patch (in pixels) to be considered valid
#' take into consideration that this is resolution dependent (default = 10000)
#' @param out_path output directory where to save the data if not returned
#' to an R variable
#' @param plot plot the processed images for reference (default = FALSE)
#' @return A nested list (data structure) listing the extracting traits together
#' with the filename of the processed file.
#' @keywords leaf traits, specific leaf area, leaf area, Gcc, Rcc
#' @export
#' @examples
#' \dontrun{
#'
#' leaf_statistics = calculate_leaf_traits(path = "/mydata")
#'
#' }

calculate_leaf_traits = function(path = "~",
                                 dpi = 300,
                                 background = "white",
                                 min_size = 10000,
                                 reference = NULL,
                                 out_path = NULL,
                                 plot = TRUE){

  # routine to calculate Gcc / Rcc
  getCover = function(img, roi, updatevalue){

    # read in the image and ROI mask
    m = (roi == updatevalue)
    m[m == 0] = NA

    # read image bands
    r = raster::raster(img, band=1)
    g = raster::raster(img, band=2)
    b = raster::raster(img, band=3)

    # mask all values which are of no
    # importance
    r_s = raster::mask(r, m, maskevalue = NA)
    g_s = raster::mask(g, m, updatevalue = NA)
    b_s = raster::mask(b, m, updatevalue = NA)
    brightness = (r_s + g_s + b_s)

    # calculate images
    gcc = g_s / brightness
    rcc = r_s / brightness

    # average across the ROI
    gcc_mean = raster::cellStats(gcc, stat = mean)
    rcc_mean = raster::cellStats(rcc, stat = mean)

    gcc_90 = raster::quantile(gcc, probs = 0.9)
    rcc_90 = raster::quantile(rcc, probs = 0.9)

    # individual channel info
    r_m = raster::cellStats(r_s, stat = mean)
    g_m = raster::cellStats(g_s, stat = mean)
    b_m = raster::cellStats(b_s, stat = mean)

    # combine raw and calculated numbers
    data = as.numeric(list(r_m,g_m,b_m,gcc_mean,gcc_90,rcc_mean,rcc_90))
    return(data)
  }

  # list all files in the directory
  # quit if no files are detected
  files = list.files(path, "*", full.names = TRUE)

  if(length(files) == 0){
    stop("No files in specified directory!")
  }

  batch_output = lapply(files,function(file){

    # feedback
    cat(sprintf("processing image: %s \n", file))

    # read in image
    im = imager::load.image(file)

    if (!is.null(reference)){
      px = imager::Xc(im) <= reference
      # blank out reference panels
      # could be done manually as well
      if (background == "white"){
        bg = imager::imfill(dim=dim(im),val=c(1,1,1))
        msk = imager::as.cimg(px)
        im = bg*msk+(1-msk)*im
      }else{
        im[px] = 0
      }
    }

    # convert to CIELAB colour space,
    # then create a data.frame with three colour channels as columns
    d = imager::sRGBtoLab(im)
    d = as.data.frame(d, wide="c")[,3:5]
    L = imager::channel(imager::RGBtoHSL(im),3)

    # convert to raster format
    L = raster::raster(t(as.matrix(L)))
    raster::extent(L) = raster::extent(0,ncol(L),0,nrow(L))

    # run k-means with 2 centers
    km = kmeans(d,2)

    # turn cluster index into an image
    seg = imager::as.cimg(km$cluster,dim=c(dim(im)[1:2],1,1))

    # convert to raster format
    seg = raster::raster(t(as.matrix(seg)))
    raster::extent(seg) = raster::extent(0,ncol(seg),0,nrow(seg))

    # calculate zonal statistics
    zstats = as.data.frame(raster::zonal(L, seg, fun = "mean", na.rm = TRUE))

    if (background == "white"){
      seg[seg == which(zstats$mean == max(zstats$mean, na.rm = TRUE))] = NA
    } else {
      seg[seg == which(zstats$mean == min(zstats$mean, na.rm = TRUE))] = NA
    }

    # clump the data
    # kick out anything smaller than 100 px
    c = raster::clump(seg, gaps = FALSE)
    raster::extent(c) = raster::extent(0,ncol(c),0,nrow(c))

    # only locations larger than min_size pixels
    for(i in which(raster::freq(c)[,2] < min_size) ){
     c[c == i] = NA
    }

    # get centroids
    clump_id = raster::getValues(c)
    xy = raster::xyFromCell(c,1:raster::ncell(c))
    df = na.omit(data.frame(xy, clump_id))
    centroids = by(df,INDICES = df$clump_id, function(xy){
      data.frame(median(xy$x[1]),
                 median(xy$y[2]))
      })
    centroids = do.call("rbind",centroids)

    # discard noise
    pixel_counts = as.data.frame(raster::freq(c, useNA = 'no'))
    results = lapply(pixel_counts$value,function(updatevalue){
                getCover(img = file,
                         roi = c,
                         updatevalue = updatevalue)
              })

    # bind things
    pixel_stats = data.frame(basename(file),
                        centroids,
                        pixel_counts$count,
                        do.call("rbind", results))
    colnames(pixel_stats) = c("filename","x","y",
                              "count","r","g","b",
                              "gcc_mean","gcc_90",
                              "rcc_mean","rcc_90")

    # convert counts to mm^2
    pixel_stats$square_cm = pixel_stats$count / (dpi / 2.54)^2

    # plot results for feedback
    if (plot){
      if(!is.null(out_path)){
        jpeg(sprintf("%s/%s",out_path,basename(file)),500,700)
      }
      raster::plotRGB(raster::brick(file))
      text(pixel_stats$x,
           pixel_stats$y,
           adj = c(0,3),
           col = "red")

      text(pixel_stats$x,
           pixel_stats$y,
           paste("Gcc: ", round(pixel_stats$gcc_90,3)),
           adj = c(0,5),
           col = "red")

      text(pixel_stats$x,
           pixel_stats$y,
           paste("cm^2: ", round(pixel_stats$square_cm,2)),
           adj = c(0,7),
           col = "red")
      if(!is.null(out_path)){
        dev.off()
      }
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
