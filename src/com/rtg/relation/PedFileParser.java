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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Properties;

import com.rtg.reference.Sex;
import com.rtg.relation.Relationship.RelationshipType;
import com.rtg.relation.Relationship.RelationshipTypeFilter;
import com.rtg.relation.Relationship.SecondInRelationshipFilter;
import com.rtg.util.StringUtils;
import com.rtg.util.diagnostic.NoTalkbackSlimException;

/**
 * Reads PLINK PED format relationships file. The format is white-space delimited
 * columns containing the following fields:
 *
 * Family ID
 * Individual ID
 * Paternal ID
 * Maternal ID
 * Sex (1=male; 2=female; other=unknown)
 * Phenotype (-9=missing, 0=missing; 1=unaffected; 2=affected)
 *
 * We ignore Family ID. There is an alternate phenotype coding, which we do not support.
 *
 */
public final class PedFileParser {

  /** Property name used to store family ID (basically informational only and for round-tripping). */
  private static final String FAMILY_ID_PROPERTY = "family-id";

  /** Reserved to indicate an unknown genome ID (only permitted to be used in maternal or paternal columns) */
  private static final String UNKNOWN = "0";

  /** Value used to denote male individuals. */
  private static final String PED_MALE = "1";

  /** Value used to denote female individuals. */
  private static final String PED_FEMALE = "2";

  /** Value used to denote "unaffected" individuals. */
  private static final String PHENOTYPE_UNAFFECTED = "1";

  /** Value used to denote "affected" individuals. */
  private static final String PHENOTYPE_AFFECTED = "2";

  /** Value used to denote individuals with unknown phenotype status. */
  private static final String PHENOTYPE_MISSING1 = "0";

  /** (Another) value used to denote individuals with unknown phenotype status. */
  private static final String PHENOTYPE_MISSING2 = "-9";

  private PedFileParser() {
  }

  /**
   * Load a PED file
   * @param file file to load
   * @return structure containing relationships and link to mapping files.
   * @throws IOException if an IO error occurs
   */
  static GenomeRelationships loadFile(File file) throws IOException {
    try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
      return load(reader);
    }
  }

  /**
   * Load a genome relationships file
   * @param reader input source
   * @return structure containing relationships and link to mapping files.
   * @throws IOException if an IO error occurs
   */
  static GenomeRelationships load(BufferedReader reader) throws IOException {

    final GenomeRelationships ped = new GenomeRelationships();
    String line;
    while ((line = reader.readLine()) != null) {
      line = line.trim(); // So we don't get thrown off by leading/trailing whitespace
      //ignore comments
      if (line.startsWith("#") || line.matches("^\\s*$")) {
        continue;
      }
      parsePedLine(ped, line);
    }

    return ped;
  }

  static String mapFromPedSex(String sex) {
    if (PED_MALE.equals(sex)) {
      return GenomeRelationships.SEX_MALE;
    } else if (PED_FEMALE.equals(sex)) {
      return GenomeRelationships.SEX_FEMALE;
    }
    return null;
  }

  static String mapToPedSex(Sex sex) {
    if (sex == Sex.MALE) {
      return PED_MALE;
    } else if (sex == Sex.FEMALE) {
      return PED_FEMALE;
    } else {
      return UNKNOWN;
    }
  }

  static String mapFromPedPhenotype(String phenotype) {
    if (PHENOTYPE_AFFECTED.equals(phenotype)) {
      return "true";
    } else if (PHENOTYPE_UNAFFECTED.equals(phenotype)) {
      return "false";
    } else if (PHENOTYPE_MISSING1.equals(phenotype)) {
      return null;
    } else if (PHENOTYPE_MISSING2.equals(phenotype)) {
      return null;
    }
    throw new NoTalkbackSlimException("Unsupported PED phenotype value: " + phenotype);
  }

  static String mapToPedPhenotype(String phenotype) {
    if ("true".equals(phenotype)) {
      return PHENOTYPE_AFFECTED;
    } else if ("false".equals(phenotype)) {
      return PHENOTYPE_UNAFFECTED;
    }
    return PHENOTYPE_MISSING1;
  }

  static void parsePedLine(GenomeRelationships ped, String line) {
    final String[] spaceSplit = line.split("\\s+");
    if (spaceSplit.length != 6) {
      throw new NoTalkbackSlimException("PED line: '" + line + "' should contain exactly 6 columns but contained: " + spaceSplit.length);
    }

    final String familyId = spaceSplit[0];
    final String individualId = spaceSplit[1];
    final String paternalId = spaceSplit[2];
    final String maternalId = spaceSplit[3];
    final String sex = spaceSplit[4];
    final String phenotype = spaceSplit[5];

    if (UNKNOWN.equals(individualId)) {
      throw new NoTalkbackSlimException("Invalid PED line: '" + line + "' Individual ID cannot be 0");
    }
    //  Specify sex if known
    ped.addGenome(individualId, mapFromPedSex(sex));

    //  Specify affliction status if known
    final String phen = mapFromPedPhenotype(phenotype);
    if (phen != null) {
      ped.getProperties(individualId).put(GenomeRelationships.DISEASE_PROPERTY, phen);
    }
    ped.getProperties(individualId).put(FAMILY_ID_PROPERTY, familyId);
    ped.getProperties(individualId).put(GenomeRelationships.PRIMARY_GENOME_PROPERTY, "true");

    if (!UNKNOWN.equals(paternalId)) {
      if (ped.getSex(paternalId) == Sex.FEMALE) { // Sex.EITHER is OK though, we'll override it now we know they are a dad
        throw new NoTalkbackSlimException("Conflicting PED definitions of sex for individual " + paternalId);
      }
      ped.addGenome(paternalId, GenomeRelationships.SEX_MALE);
      ped.addParentChild(paternalId, individualId);
    }
    if (!UNKNOWN.equals(maternalId)) {
      if (ped.getSex(maternalId) == Sex.MALE) { // Sex.EITHER is OK though, we'll override it now we know they are a mother
        throw new NoTalkbackSlimException("Conflicting PED definitions of sex for individual " + maternalId);
      }
      ped.addGenome(maternalId, GenomeRelationships.SEX_FEMALE);
      ped.addParentChild(maternalId, individualId);
    }
  }

  /**
   * Returns a PED formatted representation of the pedigree. This may not contain
   * all information in the GenomeRelationships, such as extra properties of alternate
   * relationship types.
   * @param pedigree to format as PED
   * @return the pedigree in a big old string
   */
  public static String toString(GenomeRelationships pedigree) {
    final StringBuilder sb = new StringBuilder();
    sb.append("# PED format pedigree").append(StringUtils.LS);
    sb.append("# fam-id\tind-id\tpat-id\tmat-id\tsex\tphen").append(StringUtils.LS);
    for (final String genome : pedigree.genomes()) {
      final Properties p = pedigree.getProperties(genome);
      sb.append(p.getProperty(FAMILY_ID_PROPERTY, UNKNOWN)).append("\t");
      sb.append(genome).append("\t");
      String paternal = UNKNOWN;
      String maternal = UNKNOWN;
      for (final Relationship rel : pedigree.relationships(genome, new RelationshipTypeFilter(RelationshipType.PARENT_CHILD), new SecondInRelationshipFilter(genome))) {
        final String parent = rel.first();
        if (pedigree.getSex(parent) == Sex.FEMALE) {
          maternal = parent;
        } else if (pedigree.getSex(parent) == Sex.MALE) {
          paternal = parent;
        }
      }
      sb.append(paternal).append("\t");
      sb.append(maternal).append("\t");
      sb.append(mapToPedSex(pedigree.getSex(genome))).append("\t");
      sb.append(mapToPedPhenotype(p.getProperty(GenomeRelationships.DISEASE_PROPERTY))).append(StringUtils.LS);
    }
    return sb.toString();
  }

}
