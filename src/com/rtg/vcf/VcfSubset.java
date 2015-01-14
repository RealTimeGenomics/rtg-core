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

import static com.rtg.launcher.CommonFlags.NO_GZIP;
import static com.rtg.util.cli.CommonFlagCategories.FILTERING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Set;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.vcf.header.FilterField;
import com.rtg.vcf.header.FormatField;
import com.rtg.vcf.header.IdField;
import com.rtg.vcf.header.InfoField;
import com.rtg.vcf.header.SampleField;
import com.rtg.vcf.header.VcfHeader;

/**
 */
public class VcfSubset extends AbstractCli {

  // flags
  private static final String MODULE_NAME = "vcfsubset";
  private static final String INPUT = "input";
  private static final String OUTPUT = "output";

  private static final String REMOVE_INFO = "remove-info";
  private static final String KEEP_INFO = "keep-info";
  private static final String REMOVE_INFOS = "remove-infos";

  private static final String REMOVE_FILTER = "remove-filter";
  private static final String KEEP_FILTER = "keep-filter";
  private static final String REMOVE_FILTERS = "remove-filters";

  private static final String REMOVE_SAMPLE = "remove-sample";
  private static final String KEEP_SAMPLE = "keep-sample";
  private static final String REMOVE_SAMPLES = "remove-samples";

  private static final String REMOVE_FORMAT = "remove-format";
  private static final String KEEP_FORMAT = "keep-format";

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Removes columnar data from VCF records.");
    CommonFlagCategories.setCategories(mFlags);

    mFlags.registerRequired('i', INPUT, File.class, "file", "VCF file containing variants to manipulate. Use '-' to read from standard input").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT, File.class, "file", "output VCF file. Use '-' to write to standard output").setCategory(INPUT_OUTPUT);
    CommonFlags.initNoGzip(mFlags);
    CommonFlags.initIndexFlags(mFlags);

    // Contents of FILTER
    mFlags.registerOptional(REMOVE_FILTER, String.class, "STRING", "remove the specified FILTER tag").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(KEEP_FILTER, String.class, "STRING", "keep the specified FILTER tag").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(REMOVE_FILTERS, "remove all FILTER tags").setCategory(FILTERING);

    // Contents of INFO
    mFlags.registerOptional(REMOVE_INFO, String.class, "STRING", "remove the specified INFO tag").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(KEEP_INFO, String.class, "STRING", "keep the specified INFO tag").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(REMOVE_INFOS, "remove all INFO tags").setCategory(FILTERING);

    // Contents of SAMPLE
    mFlags.registerOptional(REMOVE_SAMPLE, String.class, "STRING", "remove the specified sample").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(KEEP_SAMPLE, String.class, "STRING", "keep the specified sample").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(REMOVE_SAMPLES, "remove all samples").setCategory(FILTERING);

    // Contents of FORMAT
    mFlags.registerOptional(REMOVE_FORMAT, String.class, "STRING", "remove the specified FORMAT field").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);
    mFlags.registerOptional(KEEP_FORMAT, String.class, "STRING", "keep the specified FORMAT field").setCategory(FILTERING).setMinCount(0).setMaxCount(Integer.MAX_VALUE);

    mFlags.setValidator(new VcfSubsetValidator());
  }

  private static class VcfSubsetValidator implements Validator {
    @Override
    public boolean isValid(final CFlags flags) {
      final File input = (File) flags.getValue(INPUT);
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
      final File o = (File) flags.getValue(OUTPUT);
      if (!CommonFlags.isStdio(o)) {
        final File output = FileUtils.getZippedFileName(!flags.isSet(NO_GZIP), o);
        if (output.exists()) {
          flags.setParseMessage("The file \"" + output + "\" already exists. Please remove it first or choose a different file");
          return false;
        }
      }

      if (!checkPairNands(flags, REMOVE_INFOS, REMOVE_INFO, KEEP_INFO)) {
        return false;
      }
      if (!checkPairNands(flags, REMOVE_FILTERS, REMOVE_FILTER, KEEP_FILTER)) {
        return false;
      }
      if (!checkPairNands(flags, REMOVE_SAMPLES, REMOVE_SAMPLE, KEEP_SAMPLE)) {
        return false;
      }
      if (!flags.checkNand(REMOVE_FORMAT, KEEP_FORMAT)) {
        return false;
      }

      return true;
    }

    private boolean checkPairNands(CFlags flags, String flag1, String flag2, String flag3) {
      if (flags.isSet(flag1) && flags.isSet(flag2) || flags.isSet(flag1) && flags.isSet(flag3) || flags.isSet(flag2) && flags.isSet(flag3)) {
        flags.setParseMessage("Only one of --" + flag1 + ", --" + flag2 + ", or --" + flag3 + " can be set");
        return false;
      }
      return true;
    }
  }

  private abstract class AnnotatorAdder {

    abstract List<? extends IdField<?>> getHeaderFields(VcfHeader header);

    void additionalChecks(Set<String> flagValues, VcfHeader header) { }
    abstract VcfAnnotator makeAnnotator(boolean removeAll);
    abstract VcfAnnotator makeAnnotator(Set<String> fieldIdsSet, boolean keep);

    private VcfAnnotator processFlags(List<VcfAnnotator> annotators, VcfHeader header, String removeFlag, String keepFlag, String fieldname) {
      return processFlags(annotators, header, removeFlag, keepFlag, null, fieldname, true);
    }
    private VcfAnnotator processFlags(List<VcfAnnotator> annotators, VcfHeader header, String removeFlag, String keepFlag, String removeAllFlag, String fieldname, boolean checkHeader) {
      if (removeAllFlag != null && mFlags.isSet(removeAllFlag)) {
        final VcfAnnotator annotator = makeAnnotator(true);
        annotators.add(annotator);
        return annotator;
      } else {
        if (mFlags.isSet(removeFlag) || mFlags.isSet(keepFlag)) {
          final List<Object> infoslist;
          final boolean keep;
          if (mFlags.isSet(removeFlag)) {
            infoslist = mFlags.getValues(removeFlag);
            keep = false;
          } else {
            infoslist = mFlags.getValues(keepFlag);
            keep = true;
          }

          final Set<String> infosset = new LinkedHashSet<>();
          for (Object anInfoslist : infoslist) {
            infosset.add((String) anInfoslist);
          }

          if (checkHeader) {
            final Set<String> infoHeaderStrings = new LinkedHashSet<>();
            for (IdField<?> infoField : getHeaderFields(header)) {
              infoHeaderStrings.add(infoField.getId());
            }

            final Set<String> infossetdup = new LinkedHashSet<>(infosset);
            infossetdup.removeAll(infoHeaderStrings);

            if (infossetdup.size() > 0) {
              final StringBuilder sb = new StringBuilder();
              for (String s : infossetdup) {
                sb.append(s).append(' ');
              }
              throw new NoTalkbackSlimException(fieldname + " fields not contained in VCF meta-information: " + sb.toString().trim());
            }
          }

          additionalChecks(infosset, header);

          final VcfAnnotator annotator = makeAnnotator(infosset, keep);
          annotators.add(annotator);
          return annotator;
        }
      }
      return null;
    }
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final File input = (File) mFlags.getValue(INPUT);
    final File output = (File) mFlags.getValue(OUTPUT);
    final boolean gzip = !mFlags.isSet(NO_GZIP);
    final boolean index = !mFlags.isSet(CommonFlags.NO_INDEX);
    final boolean stdout = CommonFlags.isStdio(output);

    final List<VcfAnnotator> annotators = new ArrayList<>();

    final File vcfFile = stdout ? null : FileUtils.getZippedFileName(gzip, output);
    try (VcfReader reader = VcfReader.openVcfReader(input)) {
      final VcfHeader header = reader.getHeader();

      final AnnotatorAdder sampleAnnAdder = new AnnotatorAdder() {
        @Override
        List<SampleField> getHeaderFields(VcfHeader header) {
          return header.getSampleLines();
        }
        @Override
        VcfAnnotator makeAnnotator(boolean removeAll) {
          return new VcfSampleStripper(removeAll);
        }
        @Override
        VcfAnnotator makeAnnotator(Set<String> fieldIdsSet, boolean keep) {
          return new VcfSampleStripper(fieldIdsSet, keep);
        }
        @Override
        void additionalChecks(Set<String> flagValues, VcfHeader header) {
          boolean fail = false;
          final StringBuilder sb = new StringBuilder();
          for (String value : flagValues) {
            if (!header.getSampleNames().contains(value)) {
              fail = true;
              sb.append(value).append(' ');
            }
          }
          if (fail) {
            throw new NoTalkbackSlimException("Sample fields not contained in VCF header: " + sb.toString().trim());
          }
        }
      };
      sampleAnnAdder.processFlags(annotators, header, REMOVE_SAMPLE, KEEP_SAMPLE, REMOVE_SAMPLES, "Sample", false);

      final AnnotatorAdder infoAnnAdder = new AnnotatorAdder() {
        @Override
        List<InfoField> getHeaderFields(VcfHeader header) {
          return header.getInfoLines();
        }
        @Override
        VcfAnnotator makeAnnotator(boolean removeAll) {
          return new VcfInfoStripper(removeAll);
        }
        @Override
        VcfAnnotator makeAnnotator(Set<String> fieldIdsSet, boolean keep) {
          return new VcfInfoStripper(fieldIdsSet, keep);
        }
      };
      infoAnnAdder.processFlags(annotators, header, REMOVE_INFO, KEEP_INFO, REMOVE_INFOS, "Info", true);

      final AnnotatorAdder filterAnnAdder = new AnnotatorAdder() {
        @Override
        List<FilterField> getHeaderFields(VcfHeader header) {
          return header.getFilterLines();
        }
        @Override
        VcfAnnotator makeAnnotator(boolean removeAll) {
          return new VcfFilterStripper(removeAll);
        }
        @Override
        VcfAnnotator makeAnnotator(Set<String> fieldIdsSet, boolean keep) {
          return new VcfFilterStripper(fieldIdsSet, keep);
        }
      };
      filterAnnAdder.processFlags(annotators, header, REMOVE_FILTER, KEEP_FILTER, REMOVE_FILTERS, "Filter", true);

      final AnnotatorAdder formatAnnAdder = new AnnotatorAdder() {
        @Override
        List<FormatField> getHeaderFields(VcfHeader header) {
          return header.getFormatLines();
        }
        @Override
        VcfAnnotator makeAnnotator(boolean removeAll) {
          throw new UnsupportedOperationException("Cannot remove all formats.");
        }
        @Override
        VcfAnnotator makeAnnotator(Set<String> fieldIdsSet, boolean keep) {
          return new VcfFormatStripper(fieldIdsSet, keep);
        }
      };
      final VcfFormatStripper formatStripper = (VcfFormatStripper) formatAnnAdder.processFlags(annotators, header, REMOVE_FORMAT, KEEP_FORMAT, "Format");

      int skippedRecords = 0;
      for (final VcfAnnotator annotator : annotators) {
        annotator.updateHeader(header);
      }
      if (formatStripper != null) {
        formatStripper.updateHeader(header);
      }
      header.addRunInfo();
      try (VcfWriter writer = new VcfWriter(header, vcfFile, out, gzip, index)) {
        while (reader.hasNext()) {
          final VcfRecord rec = reader.next();
          for (final VcfAnnotator annotator : annotators) {
            annotator.annotate(rec);
          }
          if (formatStripper != null) {
            formatStripper.annotate(rec);
            if (formatStripper.keepRecord()) {
              writer.write(rec);
            } else {
              skippedRecords++;
            }
          } else {
            writer.write(rec);
          }
        }
      }
      if (skippedRecords > 0) {
        Diagnostic.warning("Records skipped due to invalid sample fields: " + skippedRecords);
      }
    }
    return 0;
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }
}
