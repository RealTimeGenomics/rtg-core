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
package com.rtg.relation;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;

import com.rtg.launcher.AbstractCli;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.io.LineWriter;
import com.rtg.vcf.header.VcfHeader;

/**
 */
public class PedFilterCli extends AbstractCli {

  private static final String KEEP_PRIMARY = "keep-primary";
  private static final String REMOVE_PARENTAGE = "remove-parentage";
  private static final String VCF_OUT = "vcf";

  @Override
  public String moduleName() {
    return "pedfilter";
  }

  @Override
  protected void initFlags() {
    mFlags.setDescription("Filter and convert a pedigree file.");
    CommonFlagCategories.setCategories(mFlags);

    mFlags.registerRequired(File.class, "FILE", "the pedigree file to process, may be PED or VCF, use '-' to read from stdin").setCategory(CommonFlagCategories.INPUT_OUTPUT);

    mFlags.registerOptional(KEEP_PRIMARY, "keep only primary individuals (those with a PED individual line / VCF sample column)").setCategory(CommonFlagCategories.FILTERING);
    mFlags.registerOptional(REMOVE_PARENTAGE, "remove all parent-child relationship information").setCategory(CommonFlagCategories.FILTERING);

    mFlags.registerOptional(VCF_OUT, "output pedigree in in the form of a VCF header rather than PED").setCategory(CommonFlagCategories.REPORTING);
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final File pedFile = (File) mFlags.getAnonymousValue(0);
    GenomeRelationships pedigree = GenomeRelationships.loadGenomeRelationships(pedFile);

    if (mFlags.isSet(KEEP_PRIMARY)) {
      pedigree = pedigree.filterByGenomes(new GenomeRelationships.PrimaryGenomeFilter(pedigree));
    }
    if (mFlags.isSet(REMOVE_PARENTAGE)) {
      pedigree = pedigree.filterByRelationships(new Relationship.NotFilter(new Relationship.RelationshipTypeFilter(Relationship.RelationshipType.PARENT_CHILD)));
    }

    try (LineWriter w = new LineWriter(new OutputStreamWriter(out))) {
      if (mFlags.isSet(VCF_OUT)) {      // Output the relationships as a VCF header
        final VcfHeader header = new VcfHeader();
        header.addCommonHeader();
        VcfPedigreeParser.addPedigreeFields(header, pedigree);
        for (String sample : pedigree.filterByGenomes(new GenomeRelationships.PrimaryGenomeFilter(pedigree)).genomes()) {
          header.addSampleName(sample);
        }
        w.write(header.toString());

      } else {             // Output the relationships in ped format
        w.write(PedFileParser.toString(pedigree));

      }
    }
    return 0;
  }


  /**
   * Command line entry point
   * @param args command line params
   */
  public static void main(String[] args) {
    new PedFilterCli().mainInit(args, System.out, System.err);
  }
}
