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
import com.rtg.protein.MapXCli;
import com.rtg.protein.MapxRename;
import com.rtg.reader.Cg2Sdf;
import com.rtg.reader.SamToFastq;
import com.rtg.reader.Sdf2Cg;
import com.rtg.reader.Sdf2Quala;
import com.rtg.reader.SdfSplitter;
import com.rtg.sam.Sam2Bam;
import com.rtg.sam.SamMergeCli;
import com.rtg.sam.SamValidatorCli;
import com.rtg.segregation.SegregationCheckerCli;
import com.rtg.segregation.SegregationVcfSearch;
import com.rtg.similarity.SimilarityCli;
import com.rtg.simulation.ReadSimEvalCli;
import com.rtg.simulation.cnv.CnvSimulatorCli;
import com.rtg.simulation.genome.GenomeSimulator;
import com.rtg.simulation.reads.CgSimCli;
import com.rtg.simulation.reads.ReadSimCli;
import com.rtg.simulation.variants.ChildSampleSimulatorCli;
import com.rtg.simulation.variants.DeNovoSampleSimulatorCli;
import com.rtg.simulation.variants.PriorPopulationVariantGeneratorCli;
import com.rtg.simulation.variants.SampleReplayerCli;
import com.rtg.simulation.variants.SampleSimulatorCli;
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
import com.rtg.variant.cnv.CnvProductCli;
import com.rtg.variant.cnv.CnvStatsCli;
import com.rtg.variant.coverage.CoverageCli;
import com.rtg.variant.sv.SvToolCli;
import com.rtg.variant.sv.UnmatedAugmenterCli;
import com.rtg.variant.sv.discord.DiscordantToolCli;
import com.rtg.vcf.SnpIntersection;
import com.rtg.visualization.Aview;
import com.rtg.zooma.ZoomaNativeBuildIndexCli;
import com.rtg.zooma.ZoomaNativeMapReadsCli;
import com.rtg.zooma.ZoomaNativeMapfReadsCli;

/**
 * Commands available in RTG Core.
 */
public final class CoreCommand {

  private CoreCommand() { }

  /** For formatting complete genomics data files for use by RtgCore */
  static final Command CG2SDF = new Command(new Cg2Sdf(), CommandCategory.FORMAT, ReleaseLevel.GA);

  /** For converting SDF formatted complete genomics reads back into CG TSV format */
  static final Command SDF2CG = new Command(new Sdf2Cg(), CommandCategory.FORMAT, ReleaseLevel.BETA);

  /** For converting RtgCore's data format into FASTA/QUALA format */
  static final Command SDF2QUALA = new Command(new Sdf2Quala(), CommandCategory.FORMAT, ReleaseLevel.ALPHA);

  /** For generating k-mer count histograms and blacklists */
  static final Command HASHDIST = new Command(new HashDistCli(), CommandCategory.UTILITY, ReleaseLevel.ALPHA);

  /** Read mapping with new and old technology mixed */
  static final Command MAP = new Command(new MapCli(), CommandCategory.MAPPING, ReleaseLevel.GA);

  /** Read mapping very constrained options, geared towards selecting reads to filter */
  static final Command MAPF = new Command(new MapFCli(), CommandCategory.MAPPING, ReleaseLevel.GA);

  /** CG reads mapping */
  static final Command CGMAP = new Command(new CgMapCli(), CommandCategory.MAPPING, ReleaseLevel.GA);

  /** Build index for zooma native mapping */
  static final Command ZOOMA_BUILD = new Command(new ZoomaNativeBuildIndexCli(), CommandCategory.MAPPING, ReleaseLevel.ALPHA);

  /** Map using zooma native code */
  static final Command ZOOMA_MAP = new Command(new ZoomaNativeMapReadsCli(), CommandCategory.MAPPING, ReleaseLevel.ALPHA);

  /** Contaminant filtering using zooma native code */
  static final Command ZOOMA_MAPF = new Command(new ZoomaNativeMapfReadsCli(), CommandCategory.MAPPING, ReleaseLevel.ALPHA);

  /** Runs Ngs and Alignment*/
  static final Command MAPX = new Command(new MapXCli(), CommandCategory.PROTEIN, ReleaseLevel.GA);

  /** Assemble reads into longer contigs */
  static final Command ASSEMBLE = new Command(new AssembleCli(), CommandCategory.ASSEMBLY, ReleaseLevel.BETA);

  /** Assemble reads into longer contigs */
  static final Command ADDPACBIO = new Command(new PacBioCli(), CommandCategory.ASSEMBLY, ReleaseLevel.BETA);

  /** Preprocess SAM files for use with the SV tool */
  static final Command SVPREP = new Command(new UnmatedAugmenterCli(), CommandCategory.VARIANT, ReleaseLevel.BETA);

  /** Run the structural variant detection tool */
  static final Command SV = new Command(new SvToolCli(), CommandCategory.VARIANT, ReleaseLevel.BETA);

  /** Run the discordant read structural variant detection tool */
  static final Command DISCORD = new Command(new DiscordantToolCli(), CommandCategory.VARIANT, ReleaseLevel.BETA);

  /** Do a coverage count from mappings */
  static final Command COVERAGE = new Command(new CoverageCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs variant calling*/
  static final Command SINGLETON = new Command(new SingletonCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs mendelian family variant calling. */
  static final Command MULTI_FAMILY = new Command(new FamilyCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs somatic variant calling. */
  static final Command MULTI_SOMATIC = new Command(new SomaticCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs somatic variant calling. */
  static final Command TUMOR_ONLY = new Command(new TumorOnlyCli(), CommandCategory.VARIANT, ReleaseLevel.ALPHA);

  /** Runs the mondo population/pedigree variant caller. */
  static final Command MULTI_POPULATION = new Command(new PopulationCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs the lineage/single parent pedigree variant caller. */
  static final Command MULTI_LINEAGE = new Command(new LineageCli(), CommandCategory.VARIANT, ReleaseLevel.BETA);

  /** Runs AVR model builder. */
  static final Command AVRBUILD = new Command(new BuilderCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs AVR model builder. */
  static final Command AVRPREDICT = new Command(new PredictCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs CNV calling */
  static final Command CNV = new Command(new CnvProductCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Runs stand alone re-calibration */
  static final Command CALIBRATE = new Command(new RecalibrateCli(), CommandCategory.VARIANT, ReleaseLevel.GA);

  /** Metagenomics species analysis */
  static final Command SPECIES = new Command(new SpeciesCli(), CommandCategory.METAGENOMICS, ReleaseLevel.GA);

  /** Metagenomics analysis of multiple heterogeneous samples */
  static final Command METASNP = new Command(new MetaSnpCli(), CommandCategory.METAGENOMICS, ReleaseLevel.ALPHA);

  /** groups similar species, also produces similarity matrices */
  static final Command SIMILARITY = new Command(new SimilarityCli(), CommandCategory.METAGENOMICS, ReleaseLevel.GA);

  /** Metagenomics composition pipeline wrapper scripting module */
  static final Command COMPOSITIONMETA = new Command(new CompositionMetaPipelineCli(), CommandCategory.PIPELINES, ReleaseLevel.BETA);

  /** Metagenomics composition pipeline wrapper scripting module */
  static final Command FUNCTIONALMETA = new Command(new FunctionalMetaPipelineCli(), CommandCategory.PIPELINES, ReleaseLevel.BETA);

  /** Metagenomics composition + functional pipeline wrapper scripting module */
  static final Command METAGENOMICS = new Command(new MetagenomicsWrapperCli(), CommandCategory.PIPELINES, ReleaseLevel.BETA);

  /** Generate simulated reads */
  static final Command GENOMESIM = new Command(new GenomeSimulator(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Generate Complete Genomics style simulated reads */
  static final Command CGSIM = new Command(new CgSimCli(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Generate simulated reads */
  static final Command READSIM = new Command(new ReadSimCli(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Evaluate read mappings and produce ROC. */
  static final Command READSIMEVAL = new Command(new ReadSimEvalCli(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Generate a VCF containing population variants for a reference */
  static final Command POPSIM = new Command(new PriorPopulationVariantGeneratorCli(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Generate a VCF containing a generated genotype for a new sample according to allele frequencies */
  static final Command SAMPLESIM = new Command(new SampleSimulatorCli(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Generate a VCF containing a generated child genotype for a new sample from two parents */
  static final Command CHILDSIM = new Command(new ChildSampleSimulatorCli(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Generate a VCF containing a derived genotype containing de novo variants */
  static final Command DENOVOSIM = new Command(new DeNovoSampleSimulatorCli(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Generate a genome SDF corresponding to a genotype contained in a VCF file */
  static final Command SAMPLEREPLAY = new Command(new SampleReplayerCli(), CommandCategory.SIMULATE, ReleaseLevel.GA);

  /** Generate a copy of a genome with a bunch of CNV mutations */
  static final Command CNVSIM = new Command(new CnvSimulatorCli(), CommandCategory.SIMULATE, ReleaseLevel.ALPHA);

  /** Evaluates CNV calling accuracy on simulated CNV data */
  static final Command CNVSIMEVAL = new Command(new CnvStatsCli(), CommandCategory.SIMULATE, ReleaseLevel.ALPHA);

  /** Classic searching with positions, also supports gaps */
  static final Command AVIEW = new Command(new Aview(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Splits an SDF file into parts */
  static final Command SDFSPLIT = new Command(new SdfSplitter(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Create a BAM file from a SAM file  */
  static final Command SAM2BAM = new Command(new Sam2Bam(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Convert a paired SAM/BAM file to FASTQ */
  static final Command SAM2FASTQ = new Command(new SamToFastq(), CommandCategory.UTILITY, ReleaseLevel.ALPHA);

  /** Rename read-id field in mapping output */
  static final Command SAMMERGE = new Command(new SamMergeCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Print statistics about a sam file */
  static final Command SAMSTATS = new Command(new SamValidatorCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Rename read-id field in mapping output */
  static final Command SAMRENAME = new Command(new SamRename(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Rename read-id field in mapping output */
  static final Command MAPXRENAME = new Command(new MapxRename(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Do checking of coverage levels from calibration information */
  static final Command CHRSTATS = new Command(new ChrStatsCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** SNP intersection class */
  static final Command SNPINTERSECT = new Command(new SnpIntersection(), CommandCategory.UTILITY, ReleaseLevel.ALPHA);

  /** NCBI Taxonomy reader */
  static final Command NCBI2TAX = new Command(new NcbiTaxDumpReaderCli(), CommandCategory.UTILITY, ReleaseLevel.BETA);

  /** Taxonomy manipulation */
  static final Command TAXFILTER = new Command(new TaxFilterCli(), CommandCategory.UTILITY, ReleaseLevel.BETA);

  /** Taxonomy verifier */
  static final Command TAXSTATS = new Command(new TaxStatsCli(), CommandCategory.UTILITY, ReleaseLevel.BETA);

  /** Search for family phasing using segregation analysis */
  static final Command PHASINGSEARCH = new Command(new SegregationVcfSearch(), CommandCategory.UTILITY, ReleaseLevel.ALPHA);

  /** Evaluate calls against segregation phasing */
  static final Command PHASINGEVAL = new Command(new SegregationCheckerCli(), CommandCategory.UTILITY, ReleaseLevel.ALPHA);

  /** Runs AVR model status. */
  static final Command AVRSTATS = new Command(new AvrStatsCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** Usage tracking server */
  static final Command USAGESERVER = new Command(new UsageServerCli(), CommandCategory.UTILITY, ReleaseLevel.GA);

  /** List modules and their license status. */
  static final LicenseCommand LICENSE = new LicenseCommand();

  /** Generate help for user. */
  static final HelpCommand HELP = new HelpCommand();

  /* This field determines the display order of the commands in the help / license output */
  private static final Command[] DISPLAY_ORDER = {
    // Formatting
    ToolsCommand.FORMAT, CG2SDF, ToolsCommand.SDF2FASTA, ToolsCommand.SDF2FASTQ, ToolsCommand.SDF2SAM, SDF2QUALA, SDF2CG,

    // Mapping
    MAP, MAPF, CGMAP,

    // Zooma native mapping
    ZOOMA_BUILD, ZOOMA_MAP, ZOOMA_MAPF,

    // Protein
    MAPX,

    // Assembly
    ASSEMBLE,
    ADDPACBIO,

    // Post-mapping
    CALIBRATE, SVPREP, SV, DISCORD, COVERAGE,

    // Variant calling
    SINGLETON, MULTI_FAMILY, MULTI_SOMATIC, MULTI_POPULATION, TUMOR_ONLY,
    MULTI_LINEAGE,
    AVRBUILD, AVRPREDICT,
    CNV,

    // Metagenomics
    SPECIES, SIMILARITY, METASNP,

    //Pipelines
    COMPOSITIONMETA,
    FUNCTIONALMETA,
    METAGENOMICS,

    // Simulation
    GENOMESIM,                                           // Reference simulation
    CGSIM, READSIM, READSIMEVAL,                         // Read simulation
    POPSIM, SAMPLESIM, CHILDSIM, DENOVOSIM, SAMPLEREPLAY, // Variant simulation
    CNVSIM, CNVSIMEVAL,                                  // Structural variant simulation

    // Utility
    ToolsCommand.BGZIP, ToolsCommand.INDEX, ToolsCommand.EXTRACT, AVIEW,                        // General purpose
    ToolsCommand.SDFSTATS, SDFSPLIT, ToolsCommand.SDFSUBSET, ToolsCommand.SDFSUBSEQ,            // SDF related
    SAM2BAM, SAMMERGE, SAMSTATS, SAMRENAME, MAPXRENAME,  // Mapping related
    CHRSTATS, SAM2FASTQ,
    ToolsCommand.MENDELIAN, PHASINGSEARCH, PHASINGEVAL,
    ToolsCommand.VCFSTATS, ToolsCommand.VCFMERGE,                       // VCF related
    ToolsCommand.VCFFILTER, ToolsCommand.VCFANNOTATE, ToolsCommand.VCFSUBSET, ToolsCommand.VCFEVAL, SNPINTERSECT,
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
