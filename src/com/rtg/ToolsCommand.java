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
package com.rtg;

import java.io.OutputStream;
import java.io.PrintStream;

import com.rtg.graph.RocPlotCli;
import com.rtg.reader.FormatCli;
import com.rtg.reader.Sdf2Fasta;
import com.rtg.reader.Sdf2Fastq;
import com.rtg.reader.SdfStatistics;
import com.rtg.reader.SdfSubseq;
import com.rtg.reader.SdfSubset;
import com.rtg.relation.PedFilterCli;
import com.rtg.relation.PedStatsCli;
import com.rtg.tabix.BgZip;
import com.rtg.tabix.ExtractCli;
import com.rtg.tabix.IndexerCli;
import com.rtg.util.License;
import com.rtg.variant.VcfStatsCli;
import com.rtg.variant.eval.VcfEvalCli;
import com.rtg.variant.util.MendeliannessChecker;
import com.rtg.vcf.VcfAnnotatorCli;
import com.rtg.vcf.VcfFilterCli;
import com.rtg.vcf.VcfMerge;
import com.rtg.vcf.VcfSubset;

/**
 * Commands available in the free tools product.
 */
public final class ToolsCommand {

  private ToolsCommand() { }

  /** For formatting data files for use by Slim */
  private static final Command FORMAT = new Command(new FormatCli(), CommandCategory.FORMAT, ReleaseLevel.GA);

  /** For converting Slim's data format into FASTA format */
  private static final Command SDF2FASTA = new Command(new Sdf2Fasta(), CommandCategory.FORMAT, ReleaseLevel.GA);

  /** For converting Slim's data format into FASTQ format */
  private static final Command SDF2FASTQ = new Command(new Sdf2Fastq(), CommandCategory.FORMAT, ReleaseLevel.GA);

  /** BGZips an input file (for use with index module) */
  private static final Command BGZIP = new Command(new BgZip(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Indexes our various output formats that have records based on reference position */
  private static final Command INDEX = new Command(new IndexerCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Extracts regions from indexed files */
  private static final Command EXTRACT = new Command(new ExtractCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Print statistics about prereads */
  private static final Command SDFSTATS = new Command(new SdfStatistics(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Creates a subset of an SDF file */
  private static final Command SDFSUBSET = new Command(new SdfSubset(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Creates a subset of an SDF file */
  private static final Command SDFSUBSEQ = new Command(new SdfSubseq(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Runs stand alone Mendelian checking */
  private static final Command MENDELIAN = new Command(new MendeliannessChecker(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** VCF stats class */
  private static final Command VCFSTATS = new Command(new VcfStatsCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** VCF merge class */
  private static final Command VCFMERGE = new Command(new VcfMerge(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** VCF subset class */
  private static final Command VCFSUBSET = new Command(new VcfSubset(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** SNP filter class */
  private static final Command VCFFILTER = new Command(new VcfFilterCli(), CommandCategory.UTILITY, ReleaseLevel.GA, License.LICENSE_KEY_PREFIX + "snpfilter");

  /** SNP filter class */
  private static final Command VCFANNOTATE = new Command(new VcfAnnotatorCli(), CommandCategory.UTILITY, ReleaseLevel.GA, License.LICENSE_KEY_PREFIX + "snpannotate");

  /** Evaluates variant calling accuracy on a given baseline variant set */
  private static final Command VCFEVAL = new Command(new VcfEvalCli(), CommandCategory.UTILITY, ReleaseLevel.GA, License.LICENSE_KEY_PREFIX + "snpsimeval");

  /** PED filter class */
  private static final Command PEDFILTER = new Command(new PedFilterCli(), CommandCategory.UTILITY, ReleaseLevel.BETA);

  /** PED stats class */
  private static final Command PEDSTATS = new Command(new PedStatsCli(), CommandCategory.UTILITY, ReleaseLevel.BETA);

  /** Roc plot tool */
  private static final Command ROCPLOT = new Command(new RocPlotCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Print version */
  private static final Command VERSION = new Command(null, "VERSION", CommandCategory.UTILITY, ReleaseLevel.GA, License.RTG_PROGRAM_KEY, null) {
    @Override
    public int mainInit(final String[] args, final OutputStream out, final PrintStream err) {
      return VersionCommand.mainInit(args, out);
    }
  };

  /** List modules and their license status. */
  private static final Command LICENSE = new Command(null, "LICENSE", CommandCategory.UTILITY, ReleaseLevel.GA, License.RTG_PROGRAM_KEY, null) {
    @Override
    public int mainInit(final String[] args, final OutputStream out, final PrintStream err) {
      return LicenseCommand.mainInit(out, INFO);
    }
  };

  /** Generate help for user. */
  static final Command HELP = new Command(null, "HELP", CommandCategory.UTILITY, ReleaseLevel.GA, License.RTG_PROGRAM_KEY, null) {
    @Override
    public int mainInit(final String[] args, final OutputStream out, final PrintStream err) {
      return HelpCommand.mainInit(args, out, err, INFO);
    }
  };

  /* This field determines the display order of the commands in the help / license output */
  private static final Command[] DISPLAY_ORDER = {
    // Formatting
    FORMAT, SDF2FASTA, SDF2FASTQ,

    // Utility
    BGZIP, INDEX, EXTRACT,                        // General purpose
    SDFSTATS, SDFSUBSET, SDFSUBSEQ,            // SDF related
    MENDELIAN, VCFSTATS, VCFMERGE,                       // VCF related
    VCFFILTER, VCFANNOTATE, VCFSUBSET, VCFEVAL,
    PEDFILTER, PEDSTATS,
    ROCPLOT,

    VERSION, LICENSE, HELP
  };

   /**
   * Provides access to list of commands
   */
  public static final CommandLookup INFO = new CommandLookup() {

    @Override
    public Command[] commands() {
      return DISPLAY_ORDER;
    }
  };
}
