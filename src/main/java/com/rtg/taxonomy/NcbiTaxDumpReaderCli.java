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

package com.rtg.taxonomy;

import static com.rtg.util.cli.CommonFlagCategories.INPUT_OUTPUT;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintStream;
import java.io.Writer;
import java.util.HashMap;
import java.util.HashSet;

import com.rtg.launcher.AbstractCli;
import com.rtg.launcher.CommonFlags;

/**
 * Reads names and nodes files from NCBI taxonomy database.  See <code>ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump_readme.txt</code>
 */
public class NcbiTaxDumpReaderCli extends AbstractCli {
  /**
   * Files in NCBI <code>taxdump</code> archive.
   */
  //private static final String CITATIONS_FILE = "citations.dmp";
  private static final String DELNODES_FILE = "delnodes.dmp";
  //private static final String DIVISION_FILE = "division.dmp";
  //private static final String GENCODE_FILE = "gencode.dmp";
  private static final String MERGED_FILE = "merged.dmp";
  private static final String NAMES_FILE = "names.dmp";
  private static final String NODES_FILE = "nodes.dmp";

  static class NcbiTaxDumpReader {

    private final HashSet<Integer> mDeletedIds = new HashSet<>();
    private final HashMap<Integer, String> mNames = new HashMap<>();

    private final Taxonomy mTaxonomy = new Taxonomy();

    private final PrintStream mErr;

    NcbiTaxDumpReader(File dir, PrintStream err) throws IOException {
      mErr = err;
      names(dir);
      nodes(dir);
      mergedIds(dir);
      deletedIds(dir);
    }

    private void deletedIds(File dir) throws IOException {
      try (BufferedReader deleteReader = new BufferedReader(new FileReader(new File(dir, DELNODES_FILE)))) {
        String line;
        while ((line = deleteReader.readLine()) != null) {
          final String[] parts = line.split("\\s*\\|\\s*");
          if (parts.length != 1) {
            mErr.println("Malformed line: " + line);
            continue;
          }
          mDeletedIds.add(Integer.valueOf(parts[0]));
        }
      }
    }

    private void mergedIds(File dir) throws IOException {
      try (BufferedReader mergedReader = new BufferedReader(new FileReader(new File(dir, MERGED_FILE)))) {
        String line;
        while ((line = mergedReader.readLine()) != null) {
          final String[] parts = line.split("\\s*\\|\\s*");
          if (parts.length != 2) {
            mErr.println("Malformed line: " + line);
            continue;
          }
          // line format is old-id then taxId
          final int oldId = Integer.parseInt(parts[0]);
          final int taxId = Integer.parseInt(parts[1]);

          mTaxonomy.addMergeNode(taxId, oldId);
        }
      }
    }

    private void names(File dir) throws IOException {
      try (BufferedReader namesReader = new BufferedReader(new FileReader(new File(dir, NAMES_FILE)))) {
        String line;
        while ((line = namesReader.readLine()) != null) {
          final String[] parts = line.split("\\s*\\|\\s*");
          if (parts.length != 4) {
            mErr.println("Malformed line: " + line);
            continue;
          }
          if ("scientific name".equals(parts[3])) {
            mNames.put(Integer.valueOf(parts[0]), new String(parts[1].toCharArray()));
            //System.out.println(parts[0] + "\t" + parts[1]);
          }
        }
      }
    }

    private void nodes(File dir) throws IOException {
      try (BufferedReader nodesReader = new BufferedReader(new FileReader(new File(dir, NODES_FILE)))) {
        String line;
        while ((line = nodesReader.readLine()) != null) {
          final String[] parts = line.split("\\s*\\|\\s*");
          if (parts.length < 12) {
            mErr.println("Malformed line [" + parts.length + "]: " + line);
            continue;
          }
          final int taxId = Integer.parseInt(parts[0]);
          final int parentId = Integer.parseInt(parts[1]);
          final String rank = new String(parts[2].toCharArray());
          final String name = mNames.get(taxId);

          mTaxonomy.addNode(taxId, parentId, name, rank);
        }
      }
    }

    public Taxonomy getTaxonomy() {
      return mTaxonomy;
    }
  }

  @Override
  public String moduleName() {
    return "ncbi2tax";
  }

  @Override
  public String description() {
    return "create RTG taxonomy from NCBI taxdump files";
  }

  @Override
  protected void initFlags() {
    mFlags.registerRequired(File.class, CommonFlags.DIR, "NCBI taxdump directory").setCategory(INPUT_OUTPUT);
    mFlags.setDescription("To convert the NCBI taxonomy into an RTG taxonomy, try:\n\n"
      + "  wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz\n"
      + "  tar zxf taxdump.tar.gz -C ncbitax\n"
      + "  rtg ncbi2tax ncbitax >rtgtax.tsv");
  }

  @Override
  protected int mainExec(OutputStream out, PrintStream err) throws IOException {

    final File dir = (File) mFlags.getAnonymousValue(0);
    err.print("Loading taxdump data from " + dir + "...");
    final NcbiTaxDumpReader ntr = new NcbiTaxDumpReader(dir, err);
    err.println("done.");

    final Taxonomy tax = ntr.getTaxonomy();

    if (!tax.isConsistent()) {
      err.println("Taxonomy is not complete.");
      err.println(tax.getInconsistencyReason());
      return 1;
    }

    err.print("Printing tree...");
    try (Writer outWriter = new OutputStreamWriter(out)) {
      tax.write(outWriter);
    }
    err.println("done.");
    return 0;
  }

}
