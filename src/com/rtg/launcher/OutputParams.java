/*
 * Copyright (c) 2014. Real Time Genomics Limited.
 *
 * Use of this source code is bound by the Real Time Genomics Limited Software Licence Agreement
 * for Academic Non-commercial Research Purposes only.
 *
 * If you did not receive a license accompanying this file, a copy must first be obtained by email
 * from support@realtimegenomics.com.  On downloading, using and/or continuing to use this source
 * code you accept the terms of that license agreement and any amendments to those terms that may
 * be made from time to time by Real Time Genomics Limited.
 */
package com.rtg.launcher;

import static com.rtg.launcher.BuildCommon.PROGRESS_FLAG;
import static com.rtg.launcher.BuildCommon.RESOURCE;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;

import com.rtg.util.ObjectParams;
import com.rtg.util.cli.CFlags;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.GzipAsynchOutputStream;

/**
 * Holds all the parameters needed for producing output to a directory.
 *
 */
public class OutputParams extends ObjectParams implements OutputDirParams {

  /**
   * Set flags.
   * @param flags to be set.
   */
  public static void initFlags(final CFlags flags) {
    flags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, "DIR", RESOURCE.getString("OUTPUT_DESC"));
    flags.registerOptional('P', PROGRESS_FLAG, RESOURCE.getString("PROGRESS_DESC"));
    CommonFlags.initNoGzip(flags);
  }

  private final File mOutputDir;

  /** True iff progress to be output.  */
  private final boolean mProgress;

  /** True iff output files to be zipped. */
  private final boolean mZip;

  /**
   * @param outDir the output directory
   * @param progress true if progress should be output
   * @param zip true if output should be zipped
   */
  public OutputParams(File outDir, boolean progress, boolean zip) {
    if (outDir == null) {
      throw new NullPointerException();
    }
    mOutputDir = outDir;
    mProgress = progress;
    mZip = zip;
    append(new Object[] {mOutputDir, mProgress, mZip});
  }

  /**
   * Get the progress flag.
   *
   * @return the progress flag. (true iff progress is to be output).
   */
  public boolean progress() {
    return mProgress;
  }

  /**
   * Get a stream to an output file in the output directory. This obeys any settings for whether compression should be applied.
   * @param name file name
   * @return the stream.
   * @throws IOException if an I/O error occurs.
   */
  public OutputStream outStream(final String name) throws IOException {
    if (!directory().exists() && !directory().mkdirs()) {
      throw new IOException("Unable to create directory \"" + directory().getPath() + "\"");
    }
    return FileUtils.createOutputStream(outFile(name), mZip, false, true);
  }


  /**
   * Return the file that will be used for output given a base name. This obeys any settings for whether compression should be applied.
   * @see #isCompressed()
   * @param name base name
   * @return the output file
   */
  public File outFile(final String name) {
    return file(mZip ? name + FileUtils.GZ_SUFFIX : name);
  }

  /**
   * @return true if output is to be compressed
   */
  public boolean isCompressed() {
    return mZip;
  }

  /**
   * @return true if output is to be compressed in blocked compressed format
   */
  public boolean isBlockCompressed() {
    return mZip && GzipAsynchOutputStream.BGZIP;
  }

  /**
   * Get the name of a child file in the output directory where all results are placed, independent of any compression settings.
   * @see #outFile(String)
   * @param child the name of the child.
   * @return the name of the file.
   */
  public File file(final String child) {
    return new File(mOutputDir, child);
  }

  /**
   * Get the output directory.
   * @return the output directory.
   */
  @Override
  public File directory() {
    return mOutputDir;
  }

  @Override
  public boolean closed() {
    return true;
  }

  @Override
  public void close() {
    // do nothing
  }

  @Override
  public String toString() {
    return "OutputParams"
        + " output directory=" + mOutputDir
        + " progress=" + mProgress
        + " zip=" + mZip;
  }

}

