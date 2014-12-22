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

package com.rtg.vcf;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;
import static com.rtg.util.cli.CommonFlagCategories.REPORTING;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.EnumSet;
import java.util.List;
import java.util.Locale;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.io.FileUtils;
import com.rtg.vcf.annotation.DerivedAnnotations;
import com.rtg.vcf.header.VcfHeader;

/**
 * Annotates a VCF file with contents of a BED file.  Annotations are appended as VCF INFO fields.
 *
 */
public final class VcfAnnotatorCli extends AbstractCli {

  private static final String MODULE_NAME = "vcfannotate";
  private static final String INPUT_FLAG = "input";
  private static final String BED_IDS_FLAG = "bed-ids";
  private static final String BED_INFO_FLAG = "bed-info";
  private static final String VCF_IDS_FLAG = "vcf-ids";
  private static final String FILL_AN_AC_FLAG = "fill-an-ac";
  private static final String X_DERIVED_ANNOTATIONS_FLAG = "Xderived-annotations";


  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Adds annotations to a VCF file.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerRequired('i', INPUT_FLAG, File.class, "file", "VCF file containing variants to annotate. Use '-' to read from standard input").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', CommonFlags.OUTPUT_FLAG, File.class, "file", "output VCF file name. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional(BED_IDS_FLAG, File.class, "file", "file in BED format containing variant ids in the name column to be added to the VCF id field").setCategory(REPORTING).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(BED_INFO_FLAG, File.class, "file", "file in BED format containing annotations in the name column to be added to the VCF info field").setCategory(REPORTING).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(VCF_IDS_FLAG, File.class, "file", "file in VCF format containing variant ids to be added to the VCF id field").setCategory(REPORTING).setMaxCount(Integer.MAX_VALUE);
    CommonFlags.initNoGzip(mFlags);
    CommonFlags.initIndexFlags(mFlags);
    final List<String> derivedRange = new ArrayList<>();
    for (final DerivedAnnotations derived : DerivedAnnotations.values()) {
      derivedRange.add(derived.toString());
    }
    mFlags.registerOptional(FILL_AN_AC_FLAG, "add or update the AN and AC info fields").setCategory(REPORTING);
    mFlags.registerOptional(X_DERIVED_ANNOTATIONS_FLAG, String.class, "STRING", "derived fields to add to VCF file").setParameterRange(derivedRange).setRangeList(true).setCategory(REPORTING);
    mFlags.setValidator(new SnpAnnotatorValidator());
  }

  private static class SnpAnnotatorValidator implements Validator {
    @Override
    public boolean isValid(CFlags flags) {
      final File input = (File) flags.getValue(INPUT_FLAG);
      if (!CommonFlags.isStdio(input)) {
        if (!input.exists()) {
          flags.setParseMessage("Given file \"" + input.getPath() + "\" does not exist.");
          return false;
        }
        if (input.isDirectory()) {
          flags.setParseMessage("Given file \"" + input.getPath() + "\" is a directory.");
          return false;
        }
      }
      if (!flags.checkNand(BED_IDS_FLAG, VCF_IDS_FLAG)) {
        return false;
      }
      if (flags.isSet(BED_INFO_FLAG) && !checkFileList(flags, BED_INFO_FLAG)) {
        return false;
      }
      if (flags.isSet(BED_IDS_FLAG) && !checkFileList(flags, BED_IDS_FLAG)) {
        return false;
      }
      if (flags.isSet(VCF_IDS_FLAG) && !checkFileList(flags, VCF_IDS_FLAG)) {
        return false;
      }
      final File o = (File) flags.getValue(CommonFlags.OUTPUT_FLAG);
      if (!CommonFlags.isStdio(o)) {
        final File output = FileUtils.getZippedFileName(!flags.isSet(CommonFlags.NO_GZIP), o);
        if (output.exists()) {
          flags.setParseMessage("The file \"" + output.getPath() + "\" already exists. Please remove it first or choose a different file");
          return false;
        }
      }
      return true;
    }

    private boolean checkFileList(CFlags flags, String flag) {
      final Collection<Object> files = flags.getValues(flag);
      for (final Object f : files) {
        final File file = (File) f;
        if (!file.exists()) {
          flags.setParseMessage("Given file \"" + file.getPath() + "\" does not exist.");
          return false;
        }
        if (file.isDirectory()) {
          flags.setParseMessage("Given file \"" + file.getPath() + "\" is a directory.");
          return false;
        }
      }
      return true;
    }
  }

  private Collection<File> getFiles(String flag) {
    final Collection<Object> inputFiles = mFlags.getValues(flag);
    final ArrayList<File> retFiles = new ArrayList<>();
    for (final Object file : inputFiles) {
      retFiles.add((File) file);
    }
    return retFiles;
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {

    final EnumSet<DerivedAnnotations> annos = EnumSet.noneOf(DerivedAnnotations.class);
    if (mFlags.isSet(X_DERIVED_ANNOTATIONS_FLAG)) {
      for (final Object anno : mFlags.getValues(X_DERIVED_ANNOTATIONS_FLAG)) {
        annos.add(DerivedAnnotations.valueOf(anno.toString().toUpperCase(Locale.getDefault())));
      }
    }
    if (mFlags.isSet(FILL_AN_AC_FLAG)) {
      annos.add(DerivedAnnotations.AN);
      annos.add(DerivedAnnotations.AC);
    }
    final List<VcfAnnotator> annotators = new ArrayList<>();
    if (mFlags.isSet(BED_INFO_FLAG)) {
      annotators.add(new BedVcfAnnotator(false, getFiles(BED_INFO_FLAG)));
    }
    if (mFlags.isSet(BED_IDS_FLAG)) {
      annotators.add(new BedVcfAnnotator(true, getFiles(BED_IDS_FLAG)));
    } else if (mFlags.isSet(VCF_IDS_FLAG)) {
      annotators.add(new VcfIdAnnotator(getFiles(VCF_IDS_FLAG)));
    }
    for (DerivedAnnotations anno : annos) {
      annotators.add(VcfUtils.getAnnotator(anno));
    }

    final File inputFile = (File) mFlags.getValue(INPUT_FLAG);
    final File output = (File) mFlags.getValue(CommonFlags.OUTPUT_FLAG);
    final boolean gzip = !mFlags.isSet(CommonFlags.NO_GZIP);
    final boolean index = !mFlags.isSet(CommonFlags.NO_INDEX);
    final boolean stdout = CommonFlags.isStdio(output);
    try (VcfReader reader = VcfReader.openVcfReader(inputFile)) {
      final VcfHeader header = reader.getHeader();
      for (final VcfAnnotator annotator : annotators) {
        annotator.updateHeader(header);
      }
      header.addRunInfo();
      final File vcfFile = stdout ? null : FileUtils.getZippedFileName(gzip, output);
      try (VcfWriter writer = new VcfWriter(header, vcfFile, out, gzip, index)) {
        while (reader.hasNext()) {
          final VcfRecord rec = reader.next();
          for (final VcfAnnotator annotator : annotators) {
            annotator.annotate(rec);
          }
          writer.write(rec);
        }
      }
    }
    return 0;
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  /**
   * main
   * @param args command line arguments
   */
  public static void main(String[] args) {
    new VcfAnnotatorCli().mainExit(args);
  }
}
