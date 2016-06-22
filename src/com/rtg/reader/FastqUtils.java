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
package com.rtg.reader;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.Arrays;
import java.util.List;

import com.rtg.mode.DnaUtils;
import com.rtg.util.io.FileUtils;

/**
 * Functions for working with FASTQ files
 */
public final class FastqUtils {

  private FastqUtils() { }


  private static final String[] EXTS = {".fastq", ".fq"};

  /**
   * @return list of extensions that we recognize for FASTQ files
   */
  public static List<String> extensions() {
    return Arrays.asList(EXTS);
  }

  /**
   * Convenience for dealing with a base file name that may have suffixes before the file extensions
   */
  public static class BaseFile {
    private final File mBaseFile;
    private final String mExtension;
    private final boolean mGzip;

    /**
     * @param baseFile file name up to point that suffixes and/or extensions are applied
     * @param extension extension that should be appended after suffix
     * @param gzip whether a gzip extension should be applied to end
     */
    public BaseFile(File baseFile, String extension, boolean gzip) {
      mBaseFile = baseFile;
      mExtension = extension;
      mGzip = gzip;
    }

    public File getBaseFile() {
      return mBaseFile;
    }

    public String getExtension() {
      return mExtension;
    }

    public boolean isGzip() {
      return mGzip;
    }

    /**
     * Apply suffix to base file and add all relevant extensions
     * @param suffix suffix for base file
     * @return the fully named file
     */
    public File suffixedFile(String suffix) {
      return new File(mBaseFile.getParentFile(), mBaseFile.getName() + suffix + mExtension + (mGzip ? FileUtils.GZ_SUFFIX : ""));
    }
  }

  /**
   * Takes a file and returns a FASTQ base file, removing any gzip extension and storing a FASTQ extension if found
   * @param file the source file
   * @param gzip whether output is intended to be gzipped
   * @return the base file
   */
  public static BaseFile baseFile(File file, boolean gzip) {
    return getBaseFile(file, gzip, EXTS);
  }

  private static BaseFile getBaseFile(File file, boolean gzip, String... exts) {
    File baseOutFile = file;
    if (FileUtils.isGzipFilename(baseOutFile)) { //remove gzip extension if present, we will put this back on if gzipping
      baseOutFile = FileUtils.removeExtension(baseOutFile);
    }
    String extension = FileUtils.getExtension(baseOutFile);
    if (Arrays.asList(exts).contains(extension)) {
      baseOutFile = FileUtils.removeExtension(baseOutFile); //remove extension if present, we will place back on with other suffixes
    } else {
      extension = exts[0];
    }
    return new BaseFile(baseOutFile, extension, gzip);
  }

  /**
   * Write a FASTQ sequence
   * @param w writer to output FASTQ to
   * @param seqName the name of the sequence
   * @param seqData the sequence data (encoded as per {@link com.rtg.mode.DNA})
   * @param qualityData the quality in raw phred values
   * @param length how long the sequence data is
   * @throws IOException if an IO error occurs
   */
  public static void writeFastqSequence(Writer w, String seqName, byte[] seqData, byte[] qualityData, int length) throws IOException {
    w.write("@");
    w.write(seqName);
    w.write("\n");
    w.write(DnaUtils.bytesToSequenceIncCG(seqData, 0, length));
    w.write("\n");
    w.write("+");
    //mAppend.append(name);
    w.write("\n");
    w.write(FastaUtils.rawToAsciiString(qualityData, 0, length));
    w.write("\n");
  }

  /**
   * Write a FASTQ sequence
   * @param w writer to output FASTQ to
   * @param seqName the name of the sequence
   * @param seqData the sequence data (encoded as per {@link com.rtg.mode.DNA})
   * @param qualityData the quality in raw phred values
   * @throws IOException if an IO error occurs
   */
  public static void writeFastqSequence(Writer w, String seqName, byte[] seqData, byte[] qualityData) throws IOException {
    writeFastqSequence(w, seqName, seqData, qualityData, seqData.length);
  }
}
