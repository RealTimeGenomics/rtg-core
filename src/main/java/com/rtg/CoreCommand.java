/*
 * Copyright (c) 2018. Real Time Genomics Limited.
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the
 *    distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
package com.rtg;

import com.rtg.assembler.AssembleCli;
import com.rtg.assembler.PacBioCli;
import com.rtg.blacklist.HashDistCli;
import com.rtg.calibrate.ChrStatsCli;
import com.rtg.calibrate.RecalibrateCli;
import com.rtg.metagenomics.CompositionMetaPipelineCli;
import com.rtg.metagenomics.FunctionalMetaPipelineCli;
import com.rtg.metagenomics.MetagenomicsWrapperCli;
import com.rtg.metagenomics.SpeciesCli;
import com.rtg.metagenomics.metasnp.MetaSnpCli;
import com.rtg.ngs.CgMapCli;
import com.rtg.ngs.MapCli;
import com.rtg.ngs.MapFCli;
import com.rtg.ngs.SamRename;
import com.rtg.protein.MapPCli;
import com.rtg.protein.MapXCli;
import com.rtg.protein.MapxRename;
import com.rtg.reader.Cg2Sdf;
import com.rtg.reader.PairedEndTrimCli;
import com.rtg.reader.SamToFastq;
import com.rtg.reader.Sdf2Cg;
import com.rtg.reader.Sdf2Quala;
import com.rtg.reader.SdfSplitter;
import com.rtg.reader.SoftClip2Fastq;
import com.rtg.sam.Sam2Bam;
import com.rtg.sam.SamMergeCli;
import com.rtg.sam.SamValidatorCli;
import com.rtg.sam.probe.DeProbeCli;
import com.rtg.segregation.SegregationCheckerCli;
import com.rtg.segregation.SegregationVcfSearch;
import com.rtg.similarity.SimilarityCli;
import com.rtg.simulation.ReadSimEvalCli;
import com.rtg.simulation.cnv.CnvSimulatorCli;
import com.rtg.taxonomy.NcbiTaxDumpReaderCli;
import com.rtg.taxonomy.TaxFilterCli;
import com.rtg.taxonomy.TaxStatsCli;
import com.rtg.usage.UsageServerCli;
import com.rtg.variant.avr.AvrStatsCli;
import com.rtg.variant.avr.BuilderCli;
import com.rtg.variant.avr.PredictCli;
import com.rtg.variant.bayes.multisample.cancer.SomaticCli;
import com.rtg.variant.bayes.multisample.cancer.TumorOnlyCli;
import com.rtg.variant.bayes.multisample.family.FamilyCli;
import com.rtg.variant.bayes.multisample.lineage.LineageCli;
import com.rtg.variant.bayes.multisample.population.PopulationCli;
import com.rtg.variant.bayes.multisample.singleton.SingletonCli;
import com.rtg.variant.cnv.CnvCli;
import com.rtg.variant.cnv.segment.CnvPonBuildCli;
import com.rtg.variant.cnv.segment.CnvSummaryCli;
import com.rtg.variant.cnv.segment.SegmentCli;
import com.rtg.variant.coverage.CoverageCli;
import com.rtg.variant.sv.FusionFilter;
import com.rtg.variant.sv.SvToolCli;
import com.rtg.variant.sv.UnmatedAugmenterCli;
import com.rtg.variant.sv.discord.DiscordantToolCli;
import com.rtg.vcf.SnpIntersection;
import com.rtg.visualization.Aview;
import com.rtg.zooma.ZoomaNativeBuildIndexCli;
import com.rtg.zooma.ZoomaNativeMapReadsCli;

/**
 * Commands available in RTG Core.
 */
public final class CoreCommand {

  private CoreCommand() { }

  /** For formatting complete genomics data files for use by RtgCore */
  static final Command CG2SDF = new LicensedCommand(new Cg2Sdf(), CommandCategory.FORMAT, ReleaseLevel.GA);

  /** For converting SDF formatted complete genomics reads back into CG TSV format */
  static final Command SDF2CG = new LicensedCommand(new Sdf2Cg(), CommandCategory.FORMAT, ReleaseLevel.BETA);

  /** Trimming of paired-end reads in FASTQ based on overlap */
  static final Command PETRIM = new LicensedCommand(new PairedEndTrimCli(), CommandCategory.FORMAT, ReleaseLevel.GA);

  /** For converting RtgCore's data format into FASTA/QUALA format */
  static final Command SDF2QUALA = new LicensedCommand(new Sdf2Quala(), CommandCategory.FORMAT, ReleaseLevel.ALPHA);

  /** For generating k-mer count histograms and blacklists */
  static final Command HASHDIST = new LicensedCommand(new HashDistCli(), CommandCategory.UTILITY, ReleaseLevel.BETA);

  /** Read mapping with new and old technology mixed */
  static final Command MAP = new LicensedCommand(new MapCli(), CommandCategory.MAPPING, ReleaseLevel.GA);

  /** Read mapping very constrained options, geared towards selecting reads to filter */
  static final Command MAPF = new LicensedCommand(new MapFCli(), CommandCategory.MAPPING, ReleaseLevel.GA);

  /** CG reads mapping */
  static final Command CGMAP = new LicensedCommand(new CgMapCli(), CommandCategory.MAPPING, ReleaseLevel.GA);

  /** Build index for zooma native mapping */
  static final Command ZOOMA_BUILD = new LicensedCommand(new ZoomaNativeBuildIndexCli(), CommandCategory.MAPPING, ReleaseLevel.ALPHA);

  /** Map using zooma native code */
  static final Command ZOOMA_MAP = new LicensedCommand(new ZoomaNativeMapReadsCli(), CommandCategory.MAPPING, ReleaseLevel.ALPHA);

  /** Do a coverage count from mappings */
  static final Command COVERAGE = new LicensedCommand(new CoverageCli(), CommandCategory.MAPPING, ReleaseLevel.GA);

  /** Runs stand alone re-calibration */
  static final Command CALIBRATE = new LicensedCommand(new RecalibrateCli(), CommandCategory.MAPPING, ReleaseLevel.GA);

  /** Protein search from translated DNA queries*/
  static final Command MAPX = new LicensedCommand(new MapXCli(), CommandCategory.PROTEIN, ReleaseLevel.GA);

  /** Protein search from untranslated protein queries*/
  static final Command MAPP = new LicensedCommand(new MapPCli(), CommandCategory.PROTEIN, ReleaseLevel.GA);

  /** Assemble reads into longer contigs */
  static final Command ASSEMBLE = new LicensedCommand(new AssembleCli(), CommandCategory.ASSEMBLY, ReleaseLevel.BETA);

  /** Assemble reads into longer contigs */
  static final Command ADDPACBIO = new LicensedCommand(new PacBioCli(), CommandCategory.ASSEMBLY, ReleaseLevel.BETA);


  /** Runs variant calling*/
  static final Command SINGLETON = new LicensedCommand(new SingletonCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs mendelian family variant calling. */
  static final Command MULTI_FAMILY = new LicensedCommand(new FamilyCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs somatic variant calling. */
  static final Command MULTI_SOMATIC = new LicensedCommand(new SomaticCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs somatic variant calling. */
  static final Command TUMOR_ONLY = new LicensedCommand(new TumorOnlyCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs the mondo population/pedigree variant caller. */
  static final Command MULTI_POPULATION = new LicensedCommand(new PopulationCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs the lineage/single parent pedigree variant caller. */
  static final Command MULTI_LINEAGE = new LicensedCommand(new LineageCli(), CommandCategory.VARIANT, ReleaseLevel.BETA);

  /** Runs AVR model builder. */
  static final Command AVRBUILD = new LicensedCommand(new BuilderCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs AVR model builder. */
  static final Command AVRPREDICT = new LicensedCommand(new PredictCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Search for family phasing using segregation analysis */
  static final Command PHASINGSEARCH = new LicensedCommand(new SegregationVcfSearch(), CommandCategory.VARIANT, ReleaseLevel.ALPHA);

  /** Evaluate calls against segregation phasing */
  static final Command PHASINGEVAL = new LicensedCommand(new SegregationCheckerCli(), CommandCategory.VARIANT, ReleaseLevel.ALPHA);

  /** Preprocess SAM files for use with the SV tool */
  static final Command SVPREP = new LicensedCommand(new UnmatedAugmenterCli(), CommandCategory.VARIANT, ReleaseLevel.BETA);

  /** Run the structural variant detection tool */
  static final Command SV = new LicensedCommand(new SvToolCli(), CommandCategory.VARIANT, ReleaseLevel.BETA);

  /** Run the discordant read structural variant detection tool */
  static final Command DISCORD = new LicensedCommand(new DiscordantToolCli(), CommandCategory.VARIANT, ReleaseLevel.BETA);

  static final Command FUSIONFILTER = new LicensedCommand(new FusionFilter(), CommandCategory.VARIANT, ReleaseLevel.ALPHA);

  /** Runs CNV calling */
  static final Command CNV = new LicensedCommand(new CnvCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs CNV segmentation. */
  static final Command SEGMENT = new LicensedCommand(new SegmentCli(), CommandCategory.VARIANT, ReleaseLevel.ALPHA);

  /** Runs CNV panel building. */
  static final Command CNVPANELBUILD = new LicensedCommand(new CnvPonBuildCli(), CommandCategory.VARIANT, ReleaseLevel.ALPHA);

  /** Runs CNV region summary report. */
  static final Command CNVSUMMARY = new LicensedCommand(new CnvSummaryCli(), CommandCategory.VARIANT, ReleaseLevel.ALPHA);


  /** Metagenomics species analysis */
  static final Command SPECIES = new LicensedCommand(new SpeciesCli(), CommandCategory.METAGENOMICS, ReleaseLevel.GA);

  /** Metagenomics analysis of multiple heterogeneous samples */
  static final Command METASNP = new LicensedCommand(new MetaSnpCli(), CommandCategory.METAGENOMICS, ReleaseLevel.ALPHA);

  /** groups similar species, also produces similarity matrices */
  static final Command SIMILARITY = new LicensedCommand(new SimilarityCli(), CommandCategory.METAGENOMICS, ReleaseLevel.GA);

  /** Metagenomics composition pipeline wrapper scripting module */
  static final Command COMPOSITIONMETA = new LicensedCommand(new CompositionMetaPipelineCli(), CommandCategory.PIPELINES, ReleaseLevel.BETA);

  /** Metagenomics composition pipeline wrapper scripting module */
  static final Command FUNCTIONALMETA = new LicensedCommand(new FunctionalMetaPipelineCli(), CommandCategory.PIPELINES, ReleaseLevel.BETA);

  /** Metagenomics composition + functional pipeline wrapper scripting module */
  static final Command METAGENOMICS = new LicensedCommand(new MetagenomicsWrapperCli(), CommandCategory.PIPELINES, ReleaseLevel.BETA);

  /** Evaluate read mappings and produce ROC. */
  static final Command READSIMEVAL = new LicensedCommand(new ReadSimEvalCli(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Generate a copy of a genome with a bunch of CNV mutations */
  static final Command CNVSIM = new LicensedCommand(new CnvSimulatorCli(), CommandCategory.SIMULATE, ReleaseLevel.ALPHA);

  /** Classic searching with positions, also supports gaps */
  static final Command AVIEW = new LicensedCommand(new Aview(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Splits an SDF file into parts */
  static final Command SDFSPLIT = new LicensedCommand(new SdfSplitter(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Create a BAM file from a SAM file  */
  static final Command SAM2BAM = new LicensedCommand(new Sam2Bam(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Convert a paired SAM/BAM file to FASTQ */
  static final Command SAM2FASTQ = new LicensedCommand(new SamToFastq(), CommandCategory.UTILITY, ReleaseLevel.ALPHA);

  /** Convert long soft-clipped alignments to FASTQ */
  static final Command SOFTCLIP2FASTQ = new LicensedCommand(new SoftClip2Fastq(), CommandCategory.UTILITY, ReleaseLevel.ALPHA);

  /** Rename read-id field in mapping output */
  static final Command SAMMERGE = new LicensedCommand(new SamMergeCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Print statistics about a sam file */
  static final Command SAMSTATS = new LicensedCommand(new SamValidatorCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Rename read-id field in mapping output */
  static final Command SAMRENAME = new LicensedCommand(new SamRename(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** strip probes for mapped bam */
  static final Command SAMSTRIPPROBES = new LicensedCommand(new DeProbeCli(), CommandCategory.UTILITY, ReleaseLevel.ALPHA);

  /** Rename read-id field in mapping output */
  static final Command MAPXRENAME = new LicensedCommand(new MapxRename(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Do checking of coverage levels from calibration information */
  static final Command CHRSTATS = new LicensedCommand(new ChrStatsCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** SNP intersection class */
  static final Command SNPINTERSECT = new LicensedCommand(new SnpIntersection(), CommandCategory.UTILITY, ReleaseLevel.ALPHA);

  /** NCBI Taxonomy reader */
  static final Command NCBI2TAX = new LicensedCommand(new NcbiTaxDumpReaderCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Taxonomy manipulation */
  static final Command TAXFILTER = new LicensedCommand(new TaxFilterCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Taxonomy verifier */
  static final Command TAXSTATS = new LicensedCommand(new TaxStatsCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Runs AVR model status. */
  static final Command AVRSTATS = new LicensedCommand(new AvrStatsCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Usage tracking server */
  static final Command USAGESERVER = new LicensedCommand(new UsageServerCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** List modules and their license status. */
  static final LicenseCommand LICENSE = new LicenseCommand();

  /** Generate help for user. */
  static final HelpCommand HELP = new HelpCommand();

  /* This field determines the display order of the commands in the help / license output */
  private static final Command[] DISPLAY_ORDER = {
    // Formatting
    ToolsCommand.FORMAT, CG2SDF, ToolsCommand.SDF2FASTA, ToolsCommand.SDF2FASTQ, ToolsCommand.SDF2SAM, SDF2QUALA, SDF2CG,
    ToolsCommand.FASTQTRIM, PETRIM,


    // Mapping
    MAP, MAPF, CGMAP,

    // Zooma native mapping
    ZOOMA_BUILD, ZOOMA_MAP,

    // Post-mapping
    COVERAGE, CALIBRATE,

    // Protein
    MAPX, MAPP,

    // Assembly
    ASSEMBLE, ADDPACBIO,


    // Variant calling
    SINGLETON, MULTI_FAMILY, MULTI_SOMATIC, MULTI_POPULATION, TUMOR_ONLY,
    MULTI_LINEAGE,
    AVRBUILD, AVRPREDICT,
    SVPREP, DISCORD, SV, FUSIONFILTER,
    CNV, SEGMENT, CNVSUMMARY, CNVPANELBUILD,
    PHASINGSEARCH, PHASINGEVAL,

    // Metagenomics
    SPECIES, SIMILARITY, METASNP,

    //Pipelines
    COMPOSITIONMETA,
    FUNCTIONALMETA,
    METAGENOMICS,

    // Simulation
    ToolsCommand.GENOMESIM,                                           // Reference simulation
    ToolsCommand.CGSIM, ToolsCommand.READSIM, READSIMEVAL,                         // Read simulation
    ToolsCommand.POPSIM, ToolsCommand.SAMPLESIM, ToolsCommand.CHILDSIM, ToolsCommand.DENOVOSIM, ToolsCommand.PEDSAMPLESIM, // Variant simulation
    ToolsCommand.SAMPLEREPLAY,               // Variant simulation
    CNVSIM,                                  // Structural variant simulation

    // Utility
    ToolsCommand.BGZIP, ToolsCommand.INDEX, ToolsCommand.EXTRACT, AVIEW,                        // General purpose
    ToolsCommand.SDFSTATS, SDFSPLIT, ToolsCommand.SDFSUBSET, ToolsCommand.SDFSUBSEQ,            // SDF related
    SAM2BAM, SAMMERGE, SAMSTATS, SAMRENAME, SAMSTRIPPROBES, MAPXRENAME,  // Mapping related
    CHRSTATS, SAM2FASTQ, SOFTCLIP2FASTQ,
    ToolsCommand.MENDELIAN,
    ToolsCommand.VCFSTATS, ToolsCommand.VCFMERGE,                       // VCF related
    ToolsCommand.VCFFILTER, ToolsCommand.VCFANNOTATE, ToolsCommand.VCFSUBSET,
    ToolsCommand.VCFSPLIT,
    ToolsCommand.VCFDECOMPOSE, ToolsCommand.VCFEVAL, ToolsCommand.VCF2ROCPLOT,
    ToolsCommand.SVDECOMPOSE, ToolsCommand.BNDEVAL, ToolsCommand.CNVEVAL,
    SNPINTERSECT,
    ToolsCommand.PEDFILTER, ToolsCommand.PEDSTATS,
    AVRSTATS, ToolsCommand.ROCPLOT,
    HASHDIST,

    NCBI2TAX, TAXFILTER, TAXSTATS, // Taxonomy

    USAGESERVER,

    ToolsCommand.VERSION, LICENSE, HELP
  };

  /**
   * Provides access to list of command
   */
  public static final CommandLookup INFO = new CommandLookup() {

    @Override
    public Command[] commands() {
      return DISPLAY_ORDER;
    }
  };

  static {
    LICENSE.setInfo(INFO);
    HELP.setInfo(INFO);
  }

}
