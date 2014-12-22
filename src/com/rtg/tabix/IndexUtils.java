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
package com.rtg.tabix;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.io.FileUtils;

import net.sf.samtools.util.BlockCompressedOutputStream;

/**
 * Utility methods for indexes
 */
public final class IndexUtils {

  private IndexUtils() { }


  /**
   * Turns a file into a block compressed file. This is a hack, do not use in production.
   * @param f file to ensure block compressed.
   * @return file that is block compressed
   * @throws java.io.IOException if an IO error occurs
   */
  public static File ensureBlockCompressed(File f) throws IOException {
    File ret = f;
    if (!TabixIndexer.isBlockCompressed(f)) {
      Diagnostic.info("blockcompressing file: " + f.getPath());
      final File outFile = File.createTempFile(f.getName(), "", f.getAbsoluteFile().getParentFile());
      try (InputStream is = FileUtils.createInputStream(f, false)) {
        try (BlockCompressedOutputStream bcos = new BlockCompressedOutputStream(outFile)) {
          FileUtils.streamToStream(is, bcos, 2048);
        }
      }
      final File mvFile;
      if (!FileUtils.isGzipFilename(f)) {
        mvFile = new File(f.getParentFile(), f.getName() + FileUtils.GZ_SUFFIX);
        if (!f.delete()) {
          Diagnostic.warning("failed to remove: " + f.getPath());
        }
        ret = mvFile;
      } else {
        mvFile = f;
      }
      if (!outFile.renameTo(mvFile)) {
        Diagnostic.warning("failed to rename temporary file: " + outFile.getPath() + " to: " + mvFile.getPath());
      }
    }
    return ret;
  }

  /**
   * Turns a list of files into a list of block compressed files. This is a hack, do not use in production.
   * @param files collection of files to ensure block compressed.
   * @return list of files that is block compressed
   * @throws java.io.IOException if an IO error occurs
   */
  public static List<File> ensureBlockCompressed(Collection<File> files) throws IOException {
    final ArrayList<File> toAdd = new ArrayList<>();
    for (File f : files) {
      toAdd.add(ensureBlockCompressed(f));
    }
    return toAdd;
  }
}
