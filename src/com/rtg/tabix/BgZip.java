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

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;

import net.sf.samtools.util.BlockCompressedOutputStream;

/**
 */
public class BgZip extends AbstractCli {

  private static final String STDOUT_FLAG = "stdout";
  private static final String DECOMPRESS_FLAG = "decompress";
  private static final String FORCE_FLAG = "force";
  private static final String NO_TERMINATE_FLAG = "no-terminate";
  private static final String LEVEL_FLAG = "compression-level";

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Compress a file with block gzip.");
    CommonFlagCategories.setCategories(mFlags);

    mFlags.registerRequired(File.class, CommonFlags.FILE, "file to (de)compress, use '-' for standard input").setCategory(INPUT_OUTPUT).setMaxCount(Integer.MAX_VALUE);

    mFlags.registerOptional('c', STDOUT_FLAG, "write on standard output, keep original files unchanged. Implied when using standard input").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional('d', DECOMPRESS_FLAG, "decompress").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional('f', FORCE_FLAG, "force overwrite of output file").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional(NO_TERMINATE_FLAG, "if set, do not add the block gzip termination block").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    mFlags.registerOptional('l', LEVEL_FLAG, Integer.class, "INT", "the compression level to use, between 1 (least but fast) and 9 (highest but slow)", BlockCompressedOutputStream.getDefaultCompressionLevel()).setCategory(CommonFlagCategories.INPUT_OUTPUT);

    mFlags.setValidator(new Validator() {
      @Override
      public boolean isValid(CFlags flags) {
        if (flags.isSet(LEVEL_FLAG) && ((Integer) flags.getValue(LEVEL_FLAG) < 1) || ((Integer) flags.getValue(LEVEL_FLAG) > 9)) {
          flags.setParseMessage("Compression level must be between 1 and 9");
          return false;
        }
        if (!flags.checkNand(STDOUT_FLAG, FORCE_FLAG)) {
          return false;
        }
        if (!flags.checkNand(DECOMPRESS_FLAG, NO_TERMINATE_FLAG)) {
          return false;
        }
        if (!flags.checkNand(DECOMPRESS_FLAG, LEVEL_FLAG)) {
          return false;
        }
        return true;
      }
    });

  }

  private static String getOutputFilename(File inputFile, boolean decompress) {
    return decompress ? inputFile.getPath().substring(0, inputFile.getPath().length() - 3) : inputFile.getPath() + FileUtils.GZ_SUFFIX;
  }

  @Override
  public String moduleName() {
    return "bgzip";
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {

    for (Object o : mFlags.getAnonymousValues(0)) {
      final File f = (File) o;
      final String outputFilename = getOutputFilename(f, mFlags.isSet(DECOMPRESS_FLAG));
      final boolean stdin = CommonFlags.isStdio(f);
      final boolean stdout = mFlags.isSet(STDOUT_FLAG) || stdin;
      if (!stdin) {
        if (!f.exists()) {
          throw new NoTalkbackSlimException("The specified file, \"" + f.getPath() + "\" does not exist.");
        } else if (mFlags.isSet(DECOMPRESS_FLAG) && !FileUtils.isGzipFilename(f)) {
          throw new NoTalkbackSlimException("Input file not in GZIP format");
        } else if (!mFlags.isSet(FORCE_FLAG) && !stdout) { //if we aren't forcibly overwriting files
          final File outfile = new File(outputFilename);
          if (outfile.exists()) {
            throw new NoTalkbackSlimException("Output file \"" + outfile.getPath() + "\" already exists.");
          }
        }
      }

      final InputStream is;
      if (stdin) {
        is = System.in;
      } else {
        is = mFlags.isSet(DECOMPRESS_FLAG) ? FileUtils.createGzipInputStream(f, true) : new FileInputStream(f);
      }
      OutputStream os = null;
      try {

        if (mFlags.isSet(DECOMPRESS_FLAG)) {
          os = stdout ? out : new FileOutputStream(outputFilename);
        } else {
          final File file = new File(f.getPath() + FileUtils.GZ_SUFFIX);
          os = stdout
            ? new BlockCompressedOutputStream(out, null, (Integer) mFlags.getValue(LEVEL_FLAG), !mFlags.isSet(NO_TERMINATE_FLAG))
            : new BlockCompressedOutputStream(new FileOutputStream(file), file, (Integer) mFlags.getValue(LEVEL_FLAG), !mFlags.isSet(NO_TERMINATE_FLAG));
        }
        final byte[] buf = new byte[FileUtils.BUFFER_SIZE];
        int bytesRead;
        while ((bytesRead = is.read(buf)) != -1) {
          os.write(buf, 0, bytesRead);
        }
      } finally {
        try {
          if (os != null) {
            os.close();
          }
        } finally {
          is.close();
        }
      }
      if (!stdout) {
        if (!f.delete()) {
          throw new NoTalkbackSlimException("Could not delete " + f.getPath());
        }
      }
    }
    return 0;
  }
}
