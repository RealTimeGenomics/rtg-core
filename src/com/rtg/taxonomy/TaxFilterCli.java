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

package com.rtg.taxonomy;

import static com.rtg.util.cli.CommonFlagCategories.FILTERING;
import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.Writer;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;
import com.rtg.reader.SdfWriter;
import com.rtg.reader.SequencesIterator;
import com.rtg.reader.SequencesReader;
import com.rtg.reader.SequencesReaderFactory;
import com.rtg.util.Constants;
import com.rtg.util.cli.CFlags;
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.cli.Validator;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.LineWriter;

/**
 */
public class TaxFilterCli extends AbstractCli {

  /** module name */
  public static final String MODULE_NAME = "taxfilter";

  private static final String SUBSET_FLAG = "subset";
  private static final String REMOVE_FLAG = "remove";
  private static final String INPUT_FLAG = "input";
  private static final String OUTPUT_FLAG = "output";
  private static final String RENAME_NORANK_FLAG = "rename-norank";
  private static final String PRUNE_INTERNAL_SEQUENCES_FLAG = "prune-internal-sequences";
  private static final String PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG = "prune-below-internal-sequences";

  static class TaxonomyFlagValidator implements Validator {
    /**
     * Checks if flags are good
     * @param flags the flags
     * @return true if good
     */
    @Override
    public boolean isValid(final CFlags flags) {
      if (!flags.checkOr(SUBSET_FLAG, REMOVE_FLAG, RENAME_NORANK_FLAG, PRUNE_INTERNAL_SEQUENCES_FLAG, PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG)) {
        return false;
      }
      if (!flags.checkNand(PRUNE_INTERNAL_SEQUENCES_FLAG, PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG)) {
        return false;
      }
      return true;
    }
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    Taxonomy tax = new Taxonomy();

    final File inputFile = (File) mFlags.getValue(INPUT_FLAG);
    final boolean sdf = inputFile.isDirectory();
    final File taxFile = sdf ? new File(inputFile, Taxonomy.TAXONOMY_FILE) : inputFile;

    if (!taxFile.exists()) {
      throw new NoTalkbackSlimException("Input taxonomy does not exist: " + taxFile);
    } else if (!taxFile.isFile()) {
      throw new NoTalkbackSlimException("Input taxonomy is not a file: " + taxFile);
    }

    //err.print("Reading taxonomy from " + taxFile + "...");
    try (FileInputStream fis = new FileInputStream(taxFile)) {
      tax.read(fis);
    }
    //err.println(tax.size());
    if (!tax.isConsistent()) {
      Diagnostic.error(tax.getInconsistencyReason());
      throw new NoTalkbackSlimException("Input taxonomy is not complete.");
    }

    if (mFlags.isSet(RENAME_NORANK_FLAG)) {
      renameNoRank(tax, (File) mFlags.getValue(RENAME_NORANK_FLAG));
    }

    if (mFlags.isSet(SUBSET_FLAG)) {
      tax = doSubset(tax, (File) mFlags.getValue(SUBSET_FLAG));
    }

    if (mFlags.isSet(REMOVE_FLAG)) {
      doRemoval(tax, (File) mFlags.getValue(REMOVE_FLAG));
    }

    final File outFile = (File) mFlags.getValue(OUTPUT_FLAG);
    if (sdf) {
      final File inMappingFile = new File(inputFile, SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE);
      if (!inMappingFile.isFile()) {
        throw new NoTalkbackSlimException("SDF has no " + SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE + " file");
      }
      final Map<String, Integer> sequenceLookupMap = SequenceToTaxonIds.sequenceToIds(inMappingFile);

      // Create a filtered sequence to id lookup file
      final Map<String, Integer> newLookup;
      if (mFlags.isSet(PRUNE_INTERNAL_SEQUENCES_FLAG)) {
        newLookup = pruneInternalSequences(tax, sequenceLookupMap);
      } else if (mFlags.isSet(PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG)) {
        newLookup = pruneBelowInternalSequences(tax, sequenceLookupMap);
      } else {
        newLookup = sequenceLookupMap;
      }

      writeSdf(inputFile, outFile, tax, newLookup, sequenceLookupMap);
    } else {
      try (Writer outWriter = new OutputStreamWriter(new FileOutputStream(outFile))) {
        tax.write(outWriter);
      }
    }

    return 0;
  }

  private void writeSdf(File inputFile, File outFile, Taxonomy tax, Map<String, Integer> sequenceLookup, Map<String, Integer> oldSequenceLookup) throws IOException {
    final File outMappingFile = new File(outFile, SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE);
    //err.println("Loaded taxonomy sequence lookup with " + sequenceLookupMap.keySet().size() + " sequences, " + new HashSet<Integer>(sequenceLookupMap.values()).size() + " taxon ids");
    try (SequencesReader reader = SequencesReaderFactory.createDefaultSequencesReaderCheckEmpty(inputFile)) {
      try (SdfWriter sdfWriter = new SdfWriter(outFile, Constants.MAX_FILE_SIZE, reader.getPrereadType(), reader.hasQualityData(), reader.hasNames(), reader.compressed(), reader.type())) {
        try (LineWriter lookupWriter = new LineWriter(new OutputStreamWriter(new FileOutputStream(outMappingFile)))) {
          sdfWriter.setPrereadArm(reader.getArm());
          // Copy any other attributes over?

          final byte[] data = new byte[(int) reader.maxLength()];
          final byte[] quality = new byte[(int) reader.maxLength()];
          final SequencesIterator it = reader.iterator();
          while (it.nextSequence()) {
            final String name = it.currentName();
            if (oldSequenceLookup.get(name) == null) {
              throw new NoTalkbackSlimException("Sequence " + name + " is present in the SDF but not in the " + SequenceToTaxonIds.TAXONOMY_TO_SEQUENCE_FILE);
            }
            final Integer taxId = sequenceLookup.get(name);
            if (taxId != null && tax.contains(taxId)) {
              final String fullname = it.currentFullName();
              final int length = it.currentLength();
              sdfWriter.startSequence(fullname);
              it.readCurrent(data);
              it.readCurrentQuality(quality);
              sdfWriter.write(data, quality, length);
              sdfWriter.endSequence();
              lookupWriter.writeln(taxId + "\t" + name);
            }
          }
        }
      }
    }
    try (Writer outWriter = new OutputStreamWriter(new FileOutputStream(new File(outFile, Taxonomy.TAXONOMY_FILE)))) {
      tax.write(outWriter);
    }
  }

  private Map<String, Integer> pruneBelowInternalSequences(Taxonomy tax, Map<String, Integer> sequenceLookupMap) {
    final Map<String, Integer> newLookup; // First axe off children of any taxon nodes that have sequence
    for (final Map.Entry<String, Integer> entry : sequenceLookupMap.entrySet()) {
      if (tax.contains(entry.getValue()) && !tax.get(entry.getValue()).isLeaf()) {
        for (final TaxonNode child : tax.get(entry.getValue()).getChildren()) {
          tax.removeSubTree(child.getId());
        }
      }
    }
    if (!tax.isConsistent()) {
      Diagnostic.error(tax.getInconsistencyReason());
      throw new NoTalkbackSlimException("Taxonomy post-prune is not complete.");
    }
    // Build new lookup from stuff that is still left
    newLookup = new HashMap<>();
    for (final Map.Entry<String, Integer> entry : sequenceLookupMap.entrySet()) {
      if (tax.contains(entry.getValue())) {
        assert tax.get(entry.getValue()).isLeaf();
        newLookup.put(entry.getKey(), entry.getValue());
      }
    }
    return newLookup;
  }

  private Map<String, Integer> pruneInternalSequences(Taxonomy tax, Map<String, Integer> sequenceLookupMap) {
    final Map<String, Integer> newLookup;
    newLookup = new HashMap<>();
    for (final Map.Entry<String, Integer> entry : sequenceLookupMap.entrySet()) {
      if (tax.contains(entry.getValue()) && tax.get(entry.getValue()).isLeaf()) {
        newLookup.put(entry.getKey(), entry.getValue());
      }
    }
    return newLookup;
  }

  private void doRemoval(Taxonomy tax, File removeFile) throws IOException {
    //err.print("Reading remove ids from " + removeFile + "...");
    final Set<Integer> ids = readIds(removeFile);
    //err.println(ids.size());

    //err.print("Removing sub-trees...");
    for (final int id : ids) {
      tax.removeSubTree(id);
    }
    //err.println(tax.size());

    if (!tax.isConsistent()) {
      Diagnostic.error(tax.getInconsistencyReason());
      throw new NoTalkbackSlimException("Taxonomy post-remove is not complete.");
    }
  }

  private Taxonomy doSubset(Taxonomy intax, File subsetFile) throws IOException {
    //err.print("Reading subset ids from " + subsetFile + "...");
    final Set<Integer> ids = readIds(subsetFile);
    //err.println(ids.size());

    //err.print("Extracting subset...");
    final Taxonomy tax;
    try {
      tax = intax.subset(ids);
    } catch (final IllegalArgumentException iae) {
      throw new NoTalkbackSlimException("Invalid subset list:" + iae.getMessage());
    }
    //err.println(tax.size());

    if (!tax.isConsistent()) {
      Diagnostic.error(tax.getInconsistencyReason());
      throw new NoTalkbackSlimException("Taxonomy post-subset is not complete.");
    }
    return tax;
  }

  private void renameNoRank(Taxonomy tax, File namesFile) throws IOException {
    final Map<Integer, String> ids = readIdsAndNames(namesFile);

    for (final Map.Entry<Integer, String> entry : ids.entrySet()) {
      final Integer id = entry.getKey();
      final TaxonNode node = tax.get(id);
      if (node != null) {
        if ("no rank".equals(node.getRank())) {
          node.setRank(entry.getValue());
        } else {
          Diagnostic.warning("Node " + id + " rank is not \"no rank\": " + node.getRank());
        }
      } else {
        Diagnostic.warning("Node not found in taxonomy: " + id);
      }
    }

    if (!tax.isConsistent()) {
      Diagnostic.error(tax.getInconsistencyReason());
      throw new NoTalkbackSlimException("Taxonomy post-rename is not complete.");
    }
  }

  private static Integer parseTaxonId(String str) throws IOException {
    try {
      return Integer.valueOf(str);
    } catch (final NumberFormatException nfe) {
      throw new IOException("Invalid taxon ID: " + str);
    }
  }

  private static Set<Integer> readIds(File file) throws IOException {
    final HashSet<Integer> ids = new HashSet<>();
    try (BufferedReader idsReader = new BufferedReader(new FileReader(file))) {
      String line;
      while ((line = idsReader.readLine()) != null) {
        final Integer taxId = parseTaxonId(line);
        if (ids.contains(taxId)) {
          Diagnostic.warning("Duplicate taxon ID: " + taxId);
        } else {
          ids.add(taxId);
        }
      }
    }
    return ids;
  }

  private static Map<Integer, String> readIdsAndNames(File file) throws IOException {
    final HashMap<Integer, String> ids = new HashMap<>();
    try (BufferedReader idsReader = new BufferedReader(new FileReader(file))) {
      String line;
      while ((line = idsReader.readLine()) != null) {
        final String[] parts = line.split("\\s+");
        final Integer taxId = parseTaxonId(parts[0]);
        final String name = parts[1];
        if (ids.containsKey(taxId)) {
          Diagnostic.warning("Duplicate taxon ID: " + taxId);
        } else {
          ids.put(taxId, name);
        }
      }
    }
    return ids;
  }

  @Override
  protected void initFlags() {
    initFlags(mFlags);
  }

  /**
   * set up a flags object for this module
   * @param flags the flags to set up
   */
  private void initFlags(CFlags flags) {
    flags.registerExtendedHelp();
    flags.setDescription("Reference taxonomy filtering.");
    CommonFlagCategories.setCategories(mFlags);
    flags.registerRequired('i', INPUT_FLAG, File.class, CommonFlags.FILE, "taxonomy input. May be either a taxonomy TSV file or an SDF containing taxonomy information").setCategory(INPUT_OUTPUT);
    flags.registerRequired('o', OUTPUT_FLAG, File.class, CommonFlags.FILE, "filename for output TSV or SDF").setCategory(INPUT_OUTPUT);
    flags.registerOptional('s', SUBSET_FLAG, File.class, CommonFlags.FILE, "file containing ids of nodes to include in subset").setCategory(FILTERING);
    flags.registerOptional('r', REMOVE_FLAG, File.class, CommonFlags.FILE, "file containing ids of nodes to remove").setCategory(FILTERING);
    flags.registerOptional('p', PRUNE_INTERNAL_SEQUENCES_FLAG, "when filtering an SDF, exclude sequence data from non-leaf output nodes").setCategory(FILTERING);
    flags.registerOptional('P', PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG, "when filtering an SDF, remove nodes below the first containing sequence data").setCategory(FILTERING);
    flags.registerOptional(RENAME_NORANK_FLAG, File.class, CommonFlags.FILE, "assign a rank to \"no rank\" nodes from file containing id/rank pairs").setCategory(FILTERING);
    flags.setValidator(new TaxonomyFlagValidator());
  }

  @Override
  public String moduleName() {
    return MODULE_NAME;
  }

  /**
   * @param args command line arguments
   */
  public static void main(String[] args) {
    new TaxFilterCli().mainExit(args);
  }
}
