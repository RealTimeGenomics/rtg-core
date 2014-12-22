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

import com.rtg.assembler.AssembleCli;
import com.rtg.assembler.PacBioCli;
import com.rtg.calibrate.ChrStatsCli;
import com.rtg.calibrate.RecalibrateCli;
import com.rtg.graph.RocPlotCli;
import com.rtg.metagenomics.CompositionMetaPipelineCli;
import com.rtg.metagenomics.FunctionalMetaPipelineCli;
import com.rtg.metagenomics.MetagenomicsWrapperCli;
import com.rtg.metagenomics.SpeciesCli;
import com.rtg.ngs.CgMapCli;
import com.rtg.ngs.MapCli;
import com.rtg.ngs.MapFCli;
import com.rtg.ngs.SamRename;
import com.rtg.protein.MapXCli;
import com.rtg.protein.MapxRename;
import com.rtg.reader.Cg2Sdf;
import com.rtg.reader.FormatCli;
import com.rtg.reader.Sdf2Fasta;
import com.rtg.reader.Sdf2Fastq;
import com.rtg.reader.Sdf2Quala;
import com.rtg.reader.SdfSplitter;
import com.rtg.reader.SdfStatistics;
import com.rtg.reader.SdfSubseq;
import com.rtg.reader.SdfSubset;
import com.rtg.relation.PedFilterCli;
import com.rtg.relation.PedStatsCli;
import com.rtg.sam.Sam2Bam;
import com.rtg.sam.SamMergeCli;
import com.rtg.sam.SamValidatorCli;
import com.rtg.segregation.SegregationCheckerCli;
import com.rtg.segregation.SegregationVcfSearch;
import com.rtg.similarity.SimilarityCli;
import com.rtg.simulation.ReadMappingAccuracy;
import com.rtg.simulation.cnv.CnvSimulatorCli;
import com.rtg.simulation.genome.GenomeSimulator;
import com.rtg.simulation.reads.CgSimCli;
import com.rtg.simulation.reads.ReadSimCli;
import com.rtg.simulation.snpsim.GenomeMutatorCli;
import com.rtg.simulation.variants.ChildSampleSimulatorCli;
import com.rtg.simulation.variants.DeNovoSampleSimulatorCli;
import com.rtg.simulation.variants.PriorPopulationVariantGeneratorCli;
import com.rtg.simulation.variants.SampleReplayerCli;
import com.rtg.simulation.variants.SampleSimulatorCli;
import com.rtg.tabix.BgZip;
import com.rtg.tabix.ExtractCli;
import com.rtg.tabix.IndexerCli;
import com.rtg.taxonomy.NcbiTaxDumpReaderCli;
import com.rtg.taxonomy.TaxFilterCli;
import com.rtg.taxonomy.TaxStatsCli;
import com.rtg.usage.UsageServerCli;
import com.rtg.util.License;
import com.rtg.variant.VcfStatsCli;
import com.rtg.variant.avr.AvrStatsCli;
import com.rtg.variant.avr.BuilderCli;
import com.rtg.variant.avr.PredictCli;
import com.rtg.variant.bayes.multisample.MultisampleCli;
import com.rtg.variant.bayes.multisample.cancer.SomaticCli;
import com.rtg.variant.bayes.multisample.family.FamilyCli;
import com.rtg.variant.bayes.multisample.lineage.LineageCli;
import com.rtg.variant.bayes.multisample.population.PopulationCli;
import com.rtg.variant.bayes.multisample.singleton.SingletonCli;
import com.rtg.variant.cnv.CnvProductCli;
import com.rtg.variant.cnv.CnvStatsCli;
import com.rtg.variant.coverage.CoverageCli;
import com.rtg.variant.eval.VcfEvalCli;
import com.rtg.variant.sv.SvToolCli;
import com.rtg.variant.sv.UnmatedAugmenterCli;
import com.rtg.variant.sv.discord.DiscordantToolCli;
import com.rtg.variant.util.MendeliannessChecker;
import com.rtg.vcf.SnpIntersection;
import com.rtg.vcf.VcfAnnotatorCli;
import com.rtg.vcf.VcfFilterCli;
import com.rtg.vcf.VcfMerge;
import com.rtg.vcf.VcfSubset;
import com.rtg.visualization.Aview;
import com.rtg.zooma.ZoomaNativeBuildIndexCli;
import com.rtg.zooma.ZoomaNativeMapReadsCli;
import com.rtg.zooma.ZoomaNativeMapfReadsCli;

/**
 * An enum of all the modules that form part of RTG Core
 */
public enum CoreCommand {

  /** For formatting data files for use by Slim */
  FORMAT(new Command(new FormatCli(), CommandCategory.FORMAT, ReleaseLevel.GA)),

  /** For formatting complete genomics data files for use by Slim */
  CG2SDF(new Command(new Cg2Sdf(), CommandCategory.FORMAT, ReleaseLevel.GA)),

  /** For converting Slim's data format into FASTA format */
  SDF2FASTA(new Command(new Sdf2Fasta(), CommandCategory.FORMAT, ReleaseLevel.GA)),

  /** For converting Slim's data format into FASTQ format */
  SDF2FASTQ(new Command(new Sdf2Fastq(), CommandCategory.FORMAT, ReleaseLevel.GA)),

  /** For converting Slim's data format into FASTA/QUALA format */
  SDF2QUALA(new Command(new Sdf2Quala(), CommandCategory.FORMAT, ReleaseLevel.ALPHA)),

  /** Read mapping with new and old technology mixed */
  MAP(new Command(new MapCli(), CommandCategory.MAPPING, ReleaseLevel.GA)),

  /** Read mapping very constrained options, geared towards selecting reads to filter */
  MAPF(new Command(new MapFCli(), CommandCategory.MAPPING, ReleaseLevel.GA)),

  /** CG reads mapping */
  CGMAP(new Command(new CgMapCli(), CommandCategory.MAPPING, ReleaseLevel.GA)),

  /** Build index for zooma native mapping */
  ZOOMA_BUILD(new Command(new ZoomaNativeBuildIndexCli(), CommandCategory.MAPPING, ReleaseLevel.ALPHA)),

  /** Map using zooma native code */
  ZOOMA_MAP(new Command(new ZoomaNativeMapReadsCli(), CommandCategory.MAPPING, ReleaseLevel.ALPHA)),

  /** Contaminant filtering using zooma native code */
  ZOOMA_MAPF(new Command(new ZoomaNativeMapfReadsCli(), CommandCategory.MAPPING, ReleaseLevel.ALPHA)),

  /** Runs Ngs and Alignment*/
  MAPX(new Command(new MapXCli(), CommandCategory.PROTEIN, ReleaseLevel.GA)),

  /** Assemble reads into longer contigs */
  ASSEMBLE(new Command(new AssembleCli(), CommandCategory.ASSEMBLY, ReleaseLevel.BETA)),

  /** Assemble reads into longer contigs */
  ADDPACBIO(new Command(new PacBioCli(), CommandCategory.ASSEMBLY, ReleaseLevel.BETA)),

  /** Preprocess SAM files for use with the SV tool */
  SVPREP(new Command(new UnmatedAugmenterCli(), CommandCategory.VARIANT, ReleaseLevel.BETA)),

  /** Run the structural variant detection tool */
  SV(new Command(new SvToolCli(), CommandCategory.VARIANT, ReleaseLevel.BETA)),

  /** Run the discordant read structural variant detection tool */
  DISCORD(new Command(new DiscordantToolCli(), CommandCategory.VARIANT, ReleaseLevel.BETA)),

  /** Do a coverage count from mappings */
  COVERAGE(new Command(new CoverageCli(), CommandCategory.VARIANT, ReleaseLevel.GA)),

  /** Runs variant calling*/
  SINGLETON(new Command(new SingletonCli(), CommandCategory.VARIANT, ReleaseLevel.GA)),

  /** Runs variant calling on multiple genomes. */
  MULTI_VARIANT(new Command(new MultisampleCli(), CommandCategory.VARIANT, ReleaseLevel.ALPHA)),

  /** Runs mendelian family variant calling. */
  MULTI_FAMILY(new Command(new FamilyCli(), CommandCategory.VARIANT, ReleaseLevel.GA)),

  /** Runs somatic variant calling. */
  MULTI_SOMATIC(new Command(new SomaticCli(), CommandCategory.VARIANT, ReleaseLevel.GA)),

  /** Runs the mondo population/pedigree variant caller. */
  MULTI_POPULATION(new Command(new PopulationCli(), CommandCategory.VARIANT, ReleaseLevel.GA)),

  /** Runs the lineage/single parent pedigree variant caller. */
  MULTI_LINEAGE(new Command(new LineageCli(), CommandCategory.VARIANT, ReleaseLevel.BETA)),

  /** Runs AVR model builder. */
  AVRBUILD(new Command(new BuilderCli(), CommandCategory.VARIANT, ReleaseLevel.GA)),

  /** Runs AVR model builder. */
  AVRPREDICT(new Command(new PredictCli(), CommandCategory.VARIANT, ReleaseLevel.GA)),

  /** Runs CNV calling */
  CNV(new Command(new CnvProductCli(), CommandCategory.VARIANT, ReleaseLevel.GA)),

  /** Runs stand alone re-calibration */
  CALIBRATE(new Command(new RecalibrateCli(), CommandCategory.VARIANT, ReleaseLevel.GA)),

  /** Metagenomics species analysis */
  SPECIES(new Command(new SpeciesCli(), CommandCategory.METAGENOMICS, ReleaseLevel.GA)),

  /** groups similar species, also produces similarity matrices */
  SIMILARITY(new Command(new SimilarityCli(), CommandCategory.METAGENOMICS, ReleaseLevel.GA)),

  /** Metagenomics composition pipeline wrapper scripting module */
  COMPOSITIONMETA(new Command(new CompositionMetaPipelineCli(), CommandCategory.PIPELINES, ReleaseLevel.BETA)),

  /** Metagenomics composition pipeline wrapper scripting module */
  FUNCTIONALMETA(new Command(new FunctionalMetaPipelineCli(), CommandCategory.PIPELINES, ReleaseLevel.BETA)),

  /** Metagenomics composition + functional pipeline wrapper scripting module */
  METAGENOMICS(new Command(new MetagenomicsWrapperCli(), CommandCategory.PIPELINES, ReleaseLevel.BETA)),

  /** Generate simulated reads */
  GENOMESIM(new Command(new GenomeSimulator(), CommandCategory.SIMULATE, ReleaseLevel.GA)),

  /** Generate Complete Genomics style simulated reads */
  CGSIM(new Command(new CgSimCli(), CommandCategory.SIMULATE, ReleaseLevel.GA)),

  /** Generate simulated reads */
  READSIM(new Command(new ReadSimCli(), CommandCategory.SIMULATE, ReleaseLevel.GA)),

  /** Evaluate read mappings and produce ROC. */
  READSIMEVAL(new Command(new ReadMappingAccuracy(), CommandCategory.SIMULATE, ReleaseLevel.BETA)),

  /**
   * Deprecated as it does not correctly allow us to simulate populations.
   * You should now use <code>popsim / samplesim / samplereplay</code>
   */
  SNPSIM(new Command(new GenomeMutatorCli(), CommandCategory.SIMULATE, ReleaseLevel.ALPHA)),

  /** Generate a VCF containing population variants for a reference */
  POPSIM(new Command(new PriorPopulationVariantGeneratorCli(), CommandCategory.SIMULATE, ReleaseLevel.BETA)),

  /** Generate a VCF containing a generated genotype for a new sample according to allele frequencies */
  SAMPLESIM(new Command(new SampleSimulatorCli(), CommandCategory.SIMULATE, ReleaseLevel.BETA)),

  /** Generate a VCF containing a generated child genotype for a new sample from two parents */
  CHILDSIM(new Command(new ChildSampleSimulatorCli(), CommandCategory.SIMULATE, ReleaseLevel.BETA)),

  /** Generate a VCF containing a derived genotype containing de novo variants */
  DENOVOSIM(new Command(new DeNovoSampleSimulatorCli(), CommandCategory.SIMULATE, ReleaseLevel.BETA)),

  /** Generate a genome SDF corresponding to a genotype contained in a VCF file */
  SAMPLEREPLAY(new Command(new SampleReplayerCli(), CommandCategory.SIMULATE, ReleaseLevel.BETA)),

  /** Evaluates variant calling accuracy on a given baseline variant set */
  VCFEVAL(new Command(new VcfEvalCli(), CommandCategory.UTILITY, ReleaseLevel.GA, License.LICENSE_KEY_PREFIX + "snpsimeval")),

  /** Generate a copy of a genome with a bunch of CNV mutations */
  CNVSIM(new Command(new CnvSimulatorCli(), CommandCategory.SIMULATE, ReleaseLevel.ALPHA)),

  /** Evaluates CNV calling accuracy on simulated CNV data */
  CNVSIMEVAL(new Command(new CnvStatsCli(), CommandCategory.SIMULATE, ReleaseLevel.ALPHA)),

  /** BGZips an input file (for use with index module) */
  BGZIP(new Command(new BgZip(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Indexes our various output formats that have records based on reference position */
  INDEX(new Command(new IndexerCli(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Extracts regions from indexed files */
  EXTRACT(new Command(new ExtractCli(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Classic searching with positions, also supports gaps */
  AVIEW(new Command(new Aview(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Print statistics about prereads */
  SDFSTATS(new Command(new SdfStatistics(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Splits an SDF file into parts */
  SDFSPLIT(new Command(new SdfSplitter(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Creates a subset of an SDF file */
  SDFSUBSET(new Command(new SdfSubset(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Creates a subset of an SDF file */
  SDFSUBSEQ(new Command(new SdfSubseq(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Create a BAM file from a SAM file  */
  SAM2BAM(new Command(new Sam2Bam(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Rename read-id field in mapping output */
  SAMMERGE(new Command(new SamMergeCli(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Print statistics about a sam file */
  SAMSTATS(new Command(new SamValidatorCli(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Rename read-id field in mapping output */
  SAMRENAME(new Command(new SamRename(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Rename read-id field in mapping output */
  MAPXRENAME(new Command(new MapxRename(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Do checking of coverage levels from calibration information */
  CHRSTATS(new Command(new ChrStatsCli(), CommandCategory.UTILITY, ReleaseLevel.BETA)),

  /** SNP filter class */
  VCFFILTER(new Command(new VcfFilterCli(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** SNP filter class */
  VCFANNOTATE(new Command(new VcfAnnotatorCli(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** SNP intersection class */
  SNPINTERSECT(new Command(new SnpIntersection(), CommandCategory.UTILITY, ReleaseLevel.ALPHA)),

  /** NCBI Taxonomy reader */
  NCBI2TAX(new Command(new NcbiTaxDumpReaderCli(), CommandCategory.UTILITY, ReleaseLevel.BETA)),

  /** Taxonomy manipulation */
  TAXFILTER(new Command(new TaxFilterCli(), CommandCategory.UTILITY, ReleaseLevel.BETA)),

  /** Taxonomy verifier */
  TAXSTATS(new Command(new TaxStatsCli(), CommandCategory.UTILITY, ReleaseLevel.BETA)),

  /** Runs stand alone Mendelian checking */
  MENDELIAN(new Command(new MendeliannessChecker(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Search for family phasing using segregation analysis */
  PHASINGSEARCH(new Command(new SegregationVcfSearch(), CommandCategory.UTILITY, ReleaseLevel.ALPHA)),

  /** Evaluate calls against segregation phasing */
  PHASINGEVAL(new Command(new SegregationCheckerCli(), CommandCategory.UTILITY, ReleaseLevel.ALPHA)),

  /** VCF stats class */
  VCFSTATS(new Command(new VcfStatsCli(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** VCF merge class */
  VCFMERGE(new Command(new VcfMerge(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** VCF subset class */
  VCFSUBSET(new Command(new VcfSubset(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** PED filter class */
  PEDFILTER(new Command(new PedFilterCli(), CommandCategory.UTILITY, ReleaseLevel.BETA)),

  /** PED stats class */
  PEDSTATS(new Command(new PedStatsCli(), CommandCategory.UTILITY, ReleaseLevel.BETA)),

  /** Runs AVR model status. */
  AVRSTATS(new Command(new AvrStatsCli(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Roc plot tool */
  ROCPLOT(new Command(new RocPlotCli(), CommandCategory.UTILITY, ReleaseLevel.GA)),

  /** Usage tracking server */
  USAGESERVER(new Command(new UsageServerCli(), CommandCategory.UTILITY, ReleaseLevel.BETA)),

  /** Print version */
  VERSION(new Command(null, "VERSION", CommandCategory.UTILITY, ReleaseLevel.GA, License.RTG_PROGRAM_KEY, null) {
    @Override
    public int mainInit(final String[] args, final OutputStream out, final PrintStream err) {
      return VersionCommand.mainInit(args, out);
    }
  }),

  /** List modules and their license status. */
  LICENSE(new Command(null, "LICENSE", CommandCategory.UTILITY, ReleaseLevel.GA, License.RTG_PROGRAM_KEY, null) {
    @Override
    public int mainInit(final String[] args, final OutputStream out, final PrintStream err) {
      return LicenseCommand.mainInit(out, INFO);
    }
  }),

  /** Generate help for user. */
  HELP(new Command(null, "HELP", CommandCategory.UTILITY, ReleaseLevel.GA, License.RTG_PROGRAM_KEY, null) {
    @Override
    public int mainInit(final String[] args, final OutputStream out, final PrintStream err) {
      return HelpCommand.mainInit(args, out, err, INFO);
    }
  });


  /* This field determines the display order of the commands in the help / license output */
  private static final CoreCommand[] DISPLAY_ORDER = {
    // Formatting
    FORMAT, CG2SDF, SDF2FASTA, SDF2FASTQ, SDF2QUALA,

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
    SINGLETON, MULTI_VARIANT, MULTI_FAMILY, MULTI_SOMATIC, MULTI_POPULATION,
    MULTI_LINEAGE,
    AVRBUILD, AVRPREDICT,
    CNV,

    // Metagenomics
    SPECIES, SIMILARITY,

    //Pipelines
    COMPOSITIONMETA,
    FUNCTIONALMETA,
    METAGENOMICS,

    // Simulation
    GENOMESIM,                                           // Reference simulation
    CGSIM, READSIM, READSIMEVAL,                         // Read simulation
    SNPSIM, POPSIM, SAMPLESIM, CHILDSIM, DENOVOSIM, SAMPLEREPLAY, // Variant simulation
    CNVSIM, CNVSIMEVAL,                                  // Structural variant simulation

    // Utility
    BGZIP, INDEX, EXTRACT, AVIEW,                        // General purpose
    SDFSTATS, SDFSPLIT, SDFSUBSET, SDFSUBSEQ,            // SDF related
    SAM2BAM, SAMMERGE, SAMSTATS, SAMRENAME, MAPXRENAME,  // Mapping related
    CHRSTATS,
    MENDELIAN, PHASINGSEARCH, PHASINGEVAL,
    VCFSTATS, VCFMERGE,                       // VCF related
    VCFFILTER, VCFANNOTATE, VCFSUBSET, VCFEVAL, SNPINTERSECT,
    PEDFILTER, PEDSTATS,
    AVRSTATS, ROCPLOT,

    NCBI2TAX, TAXFILTER, TAXSTATS, // Taxonomy

    USAGESERVER,

    VERSION, LICENSE, HELP
  };

  private final Command mModule;
  CoreCommand(final Command module) {
    mModule = module;
  }

  @Override
  public String toString() {
    return mModule.toString();
  }

  /**
   * @return the Command that can be executed
   */
  public Command module() {
    return mModule;
  }

  /**
   * Provides access to list of command
   */
  public static final CommandLookup INFO = new CommandLookup() {

    @Override
    public Command[] commands() {
      final Command[] modules = new Command[DISPLAY_ORDER.length];
      int i = 0;
      for (CoreCommand m : DISPLAY_ORDER) {
        modules[i++] = m.module();
      }
      return modules;
    }
  };
}
