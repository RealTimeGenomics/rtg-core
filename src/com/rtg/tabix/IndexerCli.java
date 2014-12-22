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
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collection;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.sam.BamIndexer;
import com.rtg.sam.SamUtils;
import com.rtg.util.Utils;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Flag;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.diagnostic.WarningType;

/**
 * Provides front end indexing various formats.
 */
public class IndexerCli extends AbstractCli {

  private static final String INPUT_FORMAT = "format";
  private static final String FORCE = "Xforce";

  /**
   * Supported formats for indexer
   */
  public static enum IndexFormat {
    /** <code>SAM</code> format */
    SAM,
    /** <code>BAM</code> format */
    BAM,
    /** <code>SV</code> format */
    SV,
    /** Coverage format */
    COVERAGETSV,
    /** BED format */
    BED,
    /** <code>VCF</code> format */
    VCF
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  private static void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Creates index files for block compressed TAB-delimited genome position files.");
    CommonFlagCategories.setCategories(flags);
    flags.setValidator(new IndexerValidator());
    flags.registerOptional(FORCE, "If the file isn't block compressed then perform this step as well").setCategory(CommonFlagCategories.UTILITY);
    final Flag inFlag = flags.registerRequired(File.class, "FILE", "block compressed files containing data to be indexed");
    inFlag.setCategory(CommonFlagCategories.INPUT_OUTPUT);
    inFlag.setMinCount(0);
    inFlag.setMaxCount(Integer.MAX_VALUE);
    inFlag.setPsuedoMinMaxRangeString(0, Integer.MAX_VALUE);
    final Flag listFlag = flags.registerOptional('I', CommonFlags.INPUT_LIST_FLAG, File.class, "FILE", "file containing a list of block compressed files (1 per line) containing genome position data").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.registerRequired('f', INPUT_FORMAT, IndexFormat.class, "FORMAT", "format of input to index").setCategory(CommonFlagCategories.INPUT_OUTPUT);
    flags.addRequiredSet(inFlag);
    flags.addRequiredSet(listFlag);
  }

  private static class IndexerValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      return CommonFlags.checkFileList(flags, CommonFlags.INPUT_LIST_FLAG, null, Integer.MAX_VALUE);
    }
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    int retCode = 0;
    final IndexFormat format = (IndexFormat) mFlags.getValue(INPUT_FORMAT);
    Collection<File> inputFiles = CommonFlags.getFileList(mFlags, CommonFlags.INPUT_LIST_FLAG, null, false);
    if (mFlags.isSet(FORCE)) {
      inputFiles = IndexUtils.ensureBlockCompressed(inputFiles);
    }
    for (final File f : inputFiles) {
      if (!TabixIndexer.isBlockCompressed(f)) {
        Diagnostic.warning("Cannot create index for " + f.getPath() + " as it is not in bgzip format.");
        retCode = 1;
        continue;
      }
      File indexFile =  new File(f.getParentFile(), f.getName() + TabixIndexer.TABIX_EXTENSION);
      final boolean indexExisted = indexFile.exists();
      try {
        switch (format) {
          case SAM:
            if (!SamUtils.looksLikeSam(f)) {
              Diagnostic.warning("File: " + f.getPath() + " does not have any headers, are you sure it is a SAM file?");
            }
            Diagnostic.info("Creating index for: " + f.getPath() + " (" + indexFile.getName() + ")");
            new TabixIndexer(f, indexFile).saveSamIndex();
            break;
          case SV:
            if (!Utils.isSvOutput(f)) {
              Diagnostic.warning("Cannot create index for " + f.getPath() + " as it is not a supported SV file.");
              retCode = 1;
              continue;
            }
            Diagnostic.info("Creating index for: " + f.getPath() + " (" + indexFile.getName() + ")");
            new TabixIndexer(f, indexFile).saveTsvIndex();
            break;
          case VCF:
            if (!Utils.isVcfFormat(f)) {
              Diagnostic.warning("Cannot create index for " + f.getPath() + " as it is not in VCF format.");
              retCode = 1;
              continue;
            }
            Diagnostic.info("Creating index for: " + f.getPath() + " (" + indexFile.getName() + ")");
            new TabixIndexer(f, indexFile).saveVcfIndex();
            break;
          case COVERAGETSV:
            if (!Utils.isCoverageOutput(f)) {
              Diagnostic.warning("Cannot create index for " + f.getPath() + " as it is not in the supported Coverage output format.");
              retCode = 1;
              continue;
            }
            Diagnostic.info("Creating index for: " + f.getPath() + " (" + indexFile.getName() + ")");
            new TabixIndexer(f, indexFile).saveBedIndex();
            break;
          case BED:
            Diagnostic.info("Creating index for: " + f.getPath() + " (" + indexFile.getName() + ")");
            new TabixIndexer(f, indexFile).saveBedIndex();
            break;
          case BAM:
            if (!SamUtils.isBAMFile(f)) {
              Diagnostic.warning("Cannot create index for " + f.getPath() + " as it is not a BAM file.");
              retCode = 1;
              continue;
            }
            //override file name for bam index
            indexFile =  new File(f.getParentFile(), f.getName() + BamIndexer.BAM_INDEX_EXTENSION);
            Diagnostic.info("Creating index for: " + f.getPath() + " (" + indexFile.getName() + ")");
            BamIndexer.saveBamIndex(f, indexFile);
            break;
          default:
            throw new RuntimeException();
        }
      } catch (IOException e) {
        if (indexFile.exists() && !indexExisted) {
          if (!indexFile.delete()) {
            Diagnostic.warning(WarningType.INFO_WARNING, "Could not delete file \"" + indexFile.getPath() + "\"");
          }
        }
        throw e;
      } catch (UnindexableDataException e) {
        if (indexFile.exists() && !indexExisted) {
          if (!indexFile.delete()) {
            Diagnostic.warning(WarningType.INFO_WARNING, "Could not delete file \"" + indexFile.getPath() + "\"");
          }
        }
        throw new NoTalkbackSlimException("Could not create index: " + e.getMessage());
      }
    }
    return retCode;
  }

  @Override
  public String moduleName() {
    return "index";
  }

  /**
   * run the program
   * @param args arguments as specified by usage
   */
  public static void main(String[] args) {
    new IndexerCli().mainExit(args);
  }
}
