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
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
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
import com.rtg.util.cli.CommonFlagCategories;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.io.FileUtils;
import com.rtg.util.io.LineWriter;

/**
 */
public class TaxFilterCli extends AbstractCli {

  private static final String SUBSET_FLAG = "subset";
  private static final String SUBTREE_FLAG = "subtree";
  private static final String REMOVE_FLAG = "remove";
  private static final String REMOVE_SEQUENCES_FLAG = "remove-sequences";
  private static final String INPUT_FLAG = "input";
  private static final String OUTPUT_FLAG = "output";
  private static final String RENAME_NORANK_FLAG = "rename-norank";
  private static final String PRUNE_INTERNAL_SEQUENCES_FLAG = "prune-internal-sequences";
  private static final String PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG = "prune-below-internal-sequences";

  @Override
  public String moduleName() {
    return "taxfilter";
  }

  @Override
  public String description() {
    return "create and manipulate RTG taxonomy files";
  }

  @Override
  protected void initFlags() {
    mFlags.registerExtendedHelp();
    mFlags.setDescription("Reference taxonomy filtering.");
    CommonFlagCategories.setCategories(mFlags);
    mFlags.registerRequired('i', INPUT_FLAG, File.class, CommonFlags.FILE, "taxonomy input. May be either a taxonomy TSV file or an SDF containing taxonomy information").setCategory(INPUT_OUTPUT);
    mFlags.registerRequired('o', OUTPUT_FLAG, File.class, CommonFlags.FILE, "filename for output TSV or SDF").setCategory(INPUT_OUTPUT);
    mFlags.registerOptional('s', SUBSET_FLAG, File.class, CommonFlags.FILE, "file containing ids of nodes to include in subset").setCategory(FILTERING);
    mFlags.registerOptional('S', SUBTREE_FLAG, File.class, CommonFlags.FILE, "file containing ids of nodes to include as subtrees in subset").setCategory(FILTERING);
    mFlags.registerOptional('r', REMOVE_FLAG, File.class, CommonFlags.FILE, "file containing ids of nodes to remove from the taxonomy (and sequence data, if any)").setCategory(FILTERING);
    mFlags.registerOptional('R', REMOVE_SEQUENCES_FLAG, File.class, CommonFlags.FILE, "file containing ids of nodes to remove sequence data from (if any)").setCategory(FILTERING);
    mFlags.registerOptional('p', PRUNE_INTERNAL_SEQUENCES_FLAG, "when filtering an SDF, exclude sequence data from non-leaf output nodes").setCategory(FILTERING);
    mFlags.registerOptional('P', PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG, "when filtering an SDF, remove nodes below the first containing sequence data").setCategory(FILTERING);
    mFlags.registerOptional(RENAME_NORANK_FLAG, File.class, CommonFlags.FILE, "assign a rank to \"no rank\" nodes from file containing id/rank pairs").setCategory(FILTERING);
    mFlags.setValidator(flags -> flags.checkOr(SUBSET_FLAG, SUBTREE_FLAG,
      REMOVE_FLAG, REMOVE_SEQUENCES_FLAG,
      RENAME_NORANK_FLAG,
      PRUNE_INTERNAL_SEQUENCES_FLAG, PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG)
      && flags.checkNand(PRUNE_INTERNAL_SEQUENCES_FLAG, PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG));
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {
    final File inputFile = (File) mFlags.getValue(INPUT_FLAG);
    final boolean sdf = inputFile.isDirectory();
    final File taxFile = sdf ? new File(inputFile, TaxonomyUtils.TAXONOMY_FILE) : inputFile;

    if (!taxFile.exists()) {
      throw new NoTalkbackSlimException("Input taxonomy does not exist: " + taxFile);
    } else if (!taxFile.isFile()) {
      throw new NoTalkbackSlimException("Input taxonomy is not a file: " + taxFile);
    }
    if (!sdf) {
      checkSdfOnlyFlags(REMOVE_SEQUENCES_FLAG, PRUNE_INTERNAL_SEQUENCES_FLAG, PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG);
    }

    //err.print("Reading taxonomy from " + taxFile + "...");
    final Taxonomy fullTax = new Taxonomy();
    try (InputStream fis = FileUtils.createInputStream(taxFile, false)) {
      fullTax.read(fis);
    }
    //err.println(tax.size());
    if (!fullTax.isConsistent()) {
      Diagnostic.error(fullTax.getInconsistencyReason());
      throw new NoTalkbackSlimException("Input taxonomy is not complete.");
    }

    if (mFlags.isSet(RENAME_NORANK_FLAG)) {
      renameNoRank(fullTax, (File) mFlags.getValue(RENAME_NORANK_FLAG));
    }

    final Taxonomy tax;
    if (mFlags.isSet(SUBSET_FLAG) || mFlags.isSet(SUBTREE_FLAG)) {
      tax = new Taxonomy();
      if (mFlags.isSet(SUBSET_FLAG)) {
        doSubset(tax, fullTax, (File) mFlags.getValue(SUBSET_FLAG));
      }
      if (mFlags.isSet(SUBTREE_FLAG)) {
        doSubTrees(tax, fullTax, (File) mFlags.getValue(SUBTREE_FLAG));
      }
    } else {
      tax = fullTax;
    }


    if (mFlags.isSet(REMOVE_FLAG)) {
      doRemoval(tax, (File) mFlags.getValue(REMOVE_FLAG));
    }

    final File outFile = (File) mFlags.getValue(OUTPUT_FLAG);
    if (sdf) {
      final File inMappingFile = new File(inputFile, TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE);
      if (!inMappingFile.isFile()) {
        throw new NoTalkbackSlimException("SDF has no " + TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE + " file");
      }
      final Map<String, Integer> sequenceLookupMap = SequenceToTaxonIds.sequenceToIds(inMappingFile);

      // Filter lookup to corresond to current taxonomy
      final Map<String, Integer> newLookup = new HashMap<>();
      for (final Map.Entry<String, Integer> entry : sequenceLookupMap.entrySet()) {
        if (tax.contains(entry.getValue())) {
          newLookup.put(entry.getKey(), entry.getValue());
        }
      }

      if (mFlags.isSet(REMOVE_SEQUENCES_FLAG)) {
        pruneSequenceData(tax, newLookup, (File) mFlags.getValue(REMOVE_SEQUENCES_FLAG));
      }
      if (mFlags.isSet(PRUNE_INTERNAL_SEQUENCES_FLAG)) {
        pruneInternalSequences(tax, newLookup);
      } else if (mFlags.isSet(PRUNE_BELOW_INTERNAL_SEQUENCES_FLAG)) {
        pruneBelowInternalSequences(tax, newLookup);
      }

      writeSdf(inputFile, outFile, tax, newLookup, sequenceLookupMap);
    } else {
      try (Writer outWriter = new OutputStreamWriter(new FileOutputStream(outFile))) {
        tax.write(outWriter);
      }
    }

    return 0;
  }

  private void checkSdfOnlyFlags(String... flags) {
    for (String flag : flags) {
      if (mFlags.isSet(flag)) {
        throw new NoTalkbackSlimException("The --" + flag + " option can only be specified when the input is a SDF");
      }
    }
  }

  private void writeSdf(File inputFile, File outFile, Taxonomy tax, Map<String, Integer> sequenceLookup, Map<String, Integer> oldSequenceLookup) throws IOException {
    final File outMappingFile = new File(outFile, TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE);
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
              throw new NoTalkbackSlimException("Sequence " + name + " is present in the SDF but not in the " + TaxonomyUtils.TAXONOMY_TO_SEQUENCE_FILE);
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
    try (Writer outWriter = new OutputStreamWriter(new FileOutputStream(new File(outFile, TaxonomyUtils.TAXONOMY_FILE)))) {
      tax.write(outWriter);
    }
  }

  private void pruneSequenceData(Taxonomy tax, Map<String, Integer> sequenceLookup, File idFile) throws IOException {
    final Set<Integer> ids = readIds(idFile);
    final Set<Integer> taxWithSeq = new HashSet<>();
    final Set<String> toRemove = new HashSet<>();

    // Collect names to remove corresponding to the tax ids
    for (final Map.Entry<String, Integer> entry : sequenceLookup.entrySet()) {
      if (ids.contains(entry.getValue())) {
        toRemove.add(entry.getKey());
      } else {
        taxWithSeq.add(entry.getValue());
      }
    }

    // Remove sequences from the lookup
    for (String seq : toRemove) {
      sequenceLookup.remove(seq);
    }

    // For each tax node, if it is a leaf, prune back to first parent either containing sequence data or with multiple children
    for (Integer id : ids) {
      TaxonNode node = tax.get(id);
      if (node != null && node.isLeaf()) {
        TaxonNode parent = node.getParent();
        while (parent != null && parent.getChildren().size() == 1 && !taxWithSeq.contains(parent.getId())) {
          node = parent;
          parent = node.getParent();
        }
        tax.removeSubTree(node.getId());
      }
    }
  }


  private void pruneBelowInternalSequences(Taxonomy tax, Map<String, Integer> sequenceLookup) {
    // Remove the subtrees from the taxonomy
    for (final Map.Entry<String, Integer> entry : sequenceLookup.entrySet()) {
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

    // Remove entries from the lookup
    final Set<String> toRemove = new HashSet<>();
    for (final Map.Entry<String, Integer> entry : sequenceLookup.entrySet()) {
      if (!tax.contains(entry.getValue())) {
        toRemove.add(entry.getKey());
      }
    }
    for (String seq : toRemove) {
      sequenceLookup.remove(seq);
    }
  }

  private void pruneInternalSequences(Taxonomy tax, Map<String, Integer> sequenceLookup) {
    final Set<String> toRemove = new HashSet<>();
    for (final Map.Entry<String, Integer> entry : sequenceLookup.entrySet()) {
      if (!tax.get(entry.getValue()).isLeaf()) {
        toRemove.add(entry.getKey());
      }
    }
    for (String seq : toRemove) {
      sequenceLookup.remove(seq);
    }
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

  private void doSubTrees(Taxonomy tax, Taxonomy source, File subtreesFile) throws IOException {
    //err.print("Reading subset ids from " + subsetFile + "...");
    final Set<Integer> ids = readIds(subtreesFile);
    //err.println(ids.size());

    //err.print("Extracting subset...");
    try {
      tax.addSubTrees(source, ids);
    } catch (final IllegalArgumentException iae) {
      throw new NoTalkbackSlimException("Invalid subset list:" + iae.getMessage());
    }
    //err.println(tax.size());

    if (!tax.isConsistent()) {
      Diagnostic.error(tax.getInconsistencyReason());
      throw new NoTalkbackSlimException("Taxonomy post-subset is not complete.");
    }
  }

  private void doSubset(Taxonomy tax, Taxonomy source, File subsetFile) throws IOException {
    //err.print("Reading subset ids from " + subsetFile + "...");
    final Set<Integer> ids = readIds(subsetFile);
    //err.println(ids.size());

    //err.print("Extracting subset...");
    try {
      tax.addPaths(source, ids);
    } catch (final IllegalArgumentException iae) {
      throw new NoTalkbackSlimException("Invalid subset list:" + iae.getMessage());
    }
    //err.println(tax.size());

    if (!tax.isConsistent()) {
      Diagnostic.error(tax.getInconsistencyReason());
      throw new NoTalkbackSlimException("Taxonomy post-subset is not complete.");
    }
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


}
