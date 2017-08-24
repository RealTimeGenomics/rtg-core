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
package com.rtg.assembler.graph.io;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.assembler.graph.Graph;
import com.rtg.assembler.graph.GraphFactory;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.implementation.ContigString;
import com.rtg.assembler.graph.implementation.PathArray;
import com.rtg.util.MD5Utils;
import com.rtg.util.StringUtils;
import com.rtg.util.gzip.WorkingGzipInputStream;
import com.rtg.util.store.StoreDirectory;
import com.rtg.util.store.StoreFile;

/**
 * Load a graph.
 *
 */
@TestClass({"com.rtg.assembler.graph.io.GraphReaderTest", "com.rtg.assembler.graph.io.CommonIOTest"})
public final class GraphReader {

  private static final String CONTIG = "contig";
  private static final String PATH = "path";
  private static final int MD5_LENGTH = 32;
  static final String MANIFEST = "manifest.txt";
  static final String HEADER_LEGACY = "header.txt";
  static final String HEADER = "header.tsv";
  private static final String VERSION = "0.0";
  static final String CONTIG_OVERLAP = "contigOverlap";

  private final StoreDirectory mGraphDirectory;
  private String mHeader = HEADER;
  private String mPathSuffix = ".tsv";

  private GraphReader(final StoreDirectory graphDir) {
    mGraphDirectory = graphDir;
  }

  private String validateDirectory() throws IOException {
    if (!mGraphDirectory.childExists(MANIFEST)) {
      //if (!new File(mGraphDirectory, MANIFEST).isFile()) {
      return "Missing " + MANIFEST;
    }
    if (!mGraphDirectory.childExists(HEADER_LEGACY) && !mGraphDirectory.childExists(HEADER)) {
      //if (!new File(mGraphDirectory, HEADER).isFile()) {
      return "Missing " + HEADER;
    }
    if (mGraphDirectory.childExists(HEADER_LEGACY)) {
      mHeader = HEADER_LEGACY;
      mPathSuffix = ".txt";
    }
    return null;
  }

  private BufferedReader reader(final String child) throws IOException {
    return new BufferedReader(new InputStreamReader(mGraphDirectory.child(child).inputStream()));
  }

  private Map<String, String> readManifest() throws IOException {
    final HashMap<String, String> manifest = new HashMap<>();
    //final BufferedReader r = new BufferedReader(new FileReader(new File(mGraphDirectory, MANIFEST)));
    try (final BufferedReader r = reader(MANIFEST)) {
      // Ignore any lines that don't look like valid md5 lines
      String line;
      while ((line = r.readLine()) != null) {
        try {
          manifest.put(line.substring(MD5_LENGTH + 1).trim(), line.substring(0, MD5_LENGTH));
        } catch (final IndexOutOfBoundsException e) {
          throw new IOException(MANIFEST + " invalid line: " + line, e);
        }
      }
    }
    return manifest;
  }

  private void checkDigest(final String file, final String md5sum, final byte[] digest) throws IOException {
    // Expect md5sum string to be hexits
    final String actual = MD5Utils.digestToString(digest);
    if (!md5sum.equals(actual)) {
      throw new IOException(file + ": md5sum error, expected " + md5sum + " got " + actual);
    }
  }

  private String nextValidLine(final BufferedReader r) throws IOException {
    while (true) {
      final String line = r.readLine();
      if (line == null) {
        return null;
      }
      final String tidy = StringUtils.removeBackslashEscapes(line);
      if (tidy.trim().length() > 0) {
        return tidy;
      }
    }
  }

  private void checkFormat(final String line, final String... token) throws IOException {
    if (line == null) {
      throw new IOException(mHeader + ": missing or incorrect date line");
    }
    final String[] parts = line.split("\t");
    if (parts.length < token.length) {
      throw new IOException(mHeader + ": expected " + token[0] + " line, saw: " + line);
    }
    for (int k = 0; k < token.length; ++k) {
      if (!parts[k].matches(token[k])) {
        throw new IOException(mHeader + ": malformed " + token[0] + " line, saw: " + line);
      }
    }
  }

  private String suckInAttributes(final Map<String, String> attr, final String attrType, final String currentLine, final BufferedReader r) throws IOException {
    String line = currentLine;
    do {
      final String[] split = line.split("\t");
      if (split.length < 2 || split.length > 3 || !split[0].equals(attrType)) {
        throw new IOException(mHeader + ": bad attribute line: " + line);
      }
      final String attrName = StringUtils.removeBackslashEscapes(split[1]);
      final String attrDescription = split.length == 3 ? StringUtils.removeBackslashEscapes(split[2]) : "";
      if (!"deleted".equals(attrName)) {
        attr.put(attrName, attrDescription);
      }
    } while ((line = nextValidLine(r)) != null && line.startsWith(attrType));
    return line;
  }
  private static final class GraphHeader {
    Map<String, String> mContigAttributes;
    Map<String, String> mPathAttributes;
    int mContigOverlap;
    private GraphHeader(int contigOverlap, Map<String, String> contigAttributes, Map<String, String> pathAttributes) {
      mContigAttributes = contigAttributes;
      mPathAttributes = pathAttributes;
      mContigOverlap = contigOverlap;
    }

  }
  private int overlap(String line) throws IOException {
    if (line == null) {
      throw new IOException(mHeader + ": missing or incorrect contigOverlap line");
    }
    final String[] parts = line.split("\t");
    if (parts.length != 2 || !CONTIG_OVERLAP.equals(parts[0])) {
      throw new IOException(mHeader + ": malformed header, expecting " + CONTIG_OVERLAP + " line, saw: " + line);
    }
    try {
      return Integer.parseInt(parts[1]);
    } catch (NumberFormatException e) {
      throw new IOException(mHeader + ": malformed header, expecting " + CONTIG_OVERLAP + " line, saw: " + line, e);
    }

  }

  // Return mappings for contig attributes and path attributes
  private GraphHeader readHeader(final String md5sum) throws IOException {
    final int contigOverlap;
    try {
      final MessageDigest md = MessageDigest.getInstance("MD5");
      final HashMap<String, String> contigAttr = new HashMap<>();
      final HashMap<String, String> pathAttr = new HashMap<>();
      try (final DigestInputStream dis = new DigestInputStream(mGraphDirectory.child(mHeader).inputStream(), md);
          final BufferedReader r = new BufferedReader(new InputStreamReader(dis))) {
        checkFormat(nextValidLine(r), "version", VERSION + ".*");
        checkFormat(nextValidLine(r), "date", "2[0-9][0-9][0-9]/[0-1][0-9]/[0-3][0-9] [0-2][0-9]:[0-5][0-9]:[0-5][0-9]");
        checkFormat(nextValidLine(r), "command", ".*");
        checkFormat(nextValidLine(r), "guid", "[0-9a-fA-F-]+");
        contigOverlap = overlap(nextValidLine(r));
        String line;
        while ((line = nextValidLine(r)) != null && line.startsWith("inputguid")) {
          checkFormat(line, "inputguid");
        }
        if (line != null) {
          line = suckInAttributes(contigAttr, CONTIG, line, r);
        }
        if (line != null) {
          suckInAttributes(pathAttr, PATH, line, r);
        }
      }
      checkDigest(mHeader, md5sum, md.digest());
      return new GraphHeader(contigOverlap, contigAttr, pathAttr);
    } catch (final NoSuchAlgorithmException e) {
      throw new RuntimeException(e);
    }
  }

  private void addContigToGraph(final MutableGraph graph, final Map<Long, Long> contigTranslationMap, final StringBuilder ntBuffer, long contigId, String attr) {
    if (contigId != 0) {
      final ContigString contig = new ContigString(ntBuffer.toString());
      final long uid = graph.addContig(contig);
      contigTranslationMap.put(contigId, uid);
      if (attr != null) {
        for (final String kv : attr.split("\t")) {
          //System.err.println("kv=" + kv);
          final int eq = kv.indexOf('=');
          final String key = kv.substring(0, eq);
          final String value = StringUtils.removeBackslashEscapes(kv.substring(eq + 1));
          if ("deleted".equals(key)) {
            if (Boolean.parseBoolean(value)) {
              graph.deleteContig(uid);
            }
          } else {
            graph.setContigAttribute(uid, key, value);
          }
        }
      }
    }
  }

  private void processContigs(final MutableGraph graph, final Map<Long, Long> contigTranslationMap, final BufferedReader r) throws IOException {
    String line;
    final StringBuilder ntBuffer = new StringBuilder();
    long contigId = 0;
    String attr = null;
    while ((line = nextValidLine(r)) != null) {
      if (line.charAt(0) == '>') {
        // Add current contig to graph before starting on the next one
        addContigToGraph(graph, contigTranslationMap, ntBuffer, contigId, attr);
        // Starting new contig
        ntBuffer.setLength(0);
        final int tab = line.indexOf('\t');
        //System.err.println(line);
        if (tab > 1) {
          contigId = Long.parseLong(line.substring(1, tab));
          attr = line.substring(tab + 1);
        } else if (tab == -1) {
          contigId = Long.parseLong(line.substring(1));
          attr = null;
        } else {
          throw new IOException("Invalid contig header line:" + line);
        }
        assert contigId > 0;
      } else {
        // Appending nucleotides to current contig
        ntBuffer.append(line);
      }
    }
    addContigToGraph(graph, contigTranslationMap, ntBuffer, contigId, attr);
  }

  private void readContigs(final MutableGraph graph, final Map<Long, Long> contigTranslationMap, final InputStream is, final boolean zipped, final String md5sum) throws IOException {
    try {
      final MessageDigest md = MessageDigest.getInstance("MD5");
      try (final DigestInputStream dis = new DigestInputStream(is, md);
          final BufferedReader r = new BufferedReader(new InputStreamReader(zipped ? new WorkingGzipInputStream(dis) : dis))) {
        processContigs(graph, contigTranslationMap, r);
      }
      checkDigest(mHeader, md5sum, md.digest());
    } catch (final NoSuchAlgorithmException e) {
      throw new RuntimeException(e);
    }
  }

  private Map<Long, Long> readContigs(final MutableGraph graph, final Map<String, String> md5sums) throws IOException {
    final Map<Long, Long> contigTranslationMap = new HashMap<>();
    int fileNumber = 0;
    while (true) {
      ++fileNumber;
      final String fname = "contig." + fileNumber + ".fa";
      if (mGraphDirectory.childExists(fname)) {
        //final File f = new File(mGraphDirectory, fname);
        //if (f.isFile()) {
        readContigs(graph, contigTranslationMap, mGraphDirectory.child(fname).inputStream(), false, md5sums.get(fname));
      } else {
        final String gname = "contig." + fileNumber + ".fa.gz";
        if (mGraphDirectory.childExists(fname)) {
          //final File g = new File(mGraphDirectory, gname);
          //if (g.isFile()) {
          readContigs(graph, contigTranslationMap, mGraphDirectory.child(gname).inputStream(), true, md5sums.get(gname));
        } else {
          // I suppose in theory should really check all files in manifest existed ...
          break;
        }
      }
    }
    return contigTranslationMap;
  }

  private void processPaths(final MutableGraph graph, final Map<Long, Long> contigTranslationMap, final BufferedReader r) throws IOException {
    String line;
    while ((line = nextValidLine(r)) != null) {
      final String[] pathSplit = line.split("\t");
      if (pathSplit.length < 2 || !PATH.equals(pathSplit[0])) {
        throw new IOException("Expected line starting with 'path', saw: " + line);
      }

      final String path = nextValidLine(r);
      if (path == null) {
        throw new IOException("Missing contigs after: " + line);
      }

      final String[] cids = path.split("\t");
      if (cids.length < 3) {
        throw new IOException("Insufficient (<2) contig identifiers on path: " + path);
      }
      if (!"contigs".equals(cids[0])) {
        throw new IOException("Expected line starting with 'contigs', saw: " + path);
      }
      final long[] translatedCids = new long[cids.length - 1];
      //System.err.println(Arrays.toString(translatedCids));
      for (int k = 0; k < translatedCids.length; ++k) {
        final String c = cids[k + 1];
        final long oldPid = Long.parseLong(c.charAt(0) == '+' ? c.substring(1) : c);
        assert oldPid != 0;
        if (oldPid > 0) {
          final Long cid = contigTranslationMap.get(oldPid);
          if (cid == null) {
            throw new IOException("Invalid contig identifier: " + oldPid + " in line:" + path);
          }
          translatedCids[k] = cid;
        } else {
          final Long cid = contigTranslationMap.get(-oldPid);
          if (cid == null) {
            throw new IOException("Invalid contig identifier: " + oldPid + " in line:" + path);
          }
          translatedCids[k] = -cid;
        }
      }
      final long pathId = graph.addPath(new PathArray(translatedCids));
      for (int i = 2; i < pathSplit.length; ++i) {
        final String kv = pathSplit[i];
        final int eq = kv.indexOf('=');
        if (eq >= 0) {
          final String key = kv.substring(0, eq);
          final String value = StringUtils.removeBackslashEscapes(kv.substring(eq + 1));
          if ("deleted".equals(key)) {
            if (Boolean.parseBoolean(value)) {
              graph.deletePath(pathId);
            }
          } else {
            graph.setPathAttribute(pathId, key, value);
          }
        }
      }
    }
  }

  private void readPaths(final MutableGraph graph, final Map<Long, Long> contigTranslationMap, final InputStream is, final boolean zipped, final String md5sum) throws IOException {
    try {
      final MessageDigest md = MessageDigest.getInstance("MD5");
      try (final DigestInputStream dis = new DigestInputStream(is, md);
          final BufferedReader r = new BufferedReader(new InputStreamReader(zipped ? new WorkingGzipInputStream(dis) : dis))) {
        processPaths(graph, contigTranslationMap, r);
      }
      checkDigest(mHeader, md5sum, md.digest());
    } catch (final NoSuchAlgorithmException e) {
      throw new RuntimeException(e);
    }
  }

  private void readPaths(final MutableGraph graph, final Map<Long, Long> contigTranslationMap, final Map<String, String> md5sums) throws IOException {
    int fileNumber = 0;
    while (true) {
      ++fileNumber;
      final String fname = "path." + fileNumber + mPathSuffix;
      if (mGraphDirectory.childExists(fname)) {
        readPaths(graph, contigTranslationMap, mGraphDirectory.child(fname).inputStream(), false, md5sums.get(fname));
      } else {
        final String gname = "path." + fileNumber + mPathSuffix + ".gz";
        if (mGraphDirectory.childExists(gname)) {
          readPaths(graph, contigTranslationMap, mGraphDirectory.child(gname).inputStream(), true, md5sums.get(gname));
        } else {
          // I suppose in theory should really check all files in manifest existed ...
          break;
        }
      }
    }
  }

  private Graph read(final GraphFactory factory, Map<String, String> contigAttributes, Map<String, String> pathAttributes) throws IOException {
    final String msg = validateDirectory();
    if (msg != null) {
      throw new IOException(msg);
    }
    final Map<String, String> md5sums = readManifest();
    final GraphHeader header = readHeader(md5sums.get(mHeader));
    if (contigAttributes != null) {
      for (Map.Entry<String, String> entry : contigAttributes.entrySet()) {
        if (!header.mContigAttributes.containsKey(entry.getKey())) {
          header.mContigAttributes.put(entry.getKey(), entry.getValue());
        }
      }
    }
    if (pathAttributes != null) {
      for (Map.Entry<String, String> entry : pathAttributes.entrySet()) {
        if (!header.mPathAttributes.containsKey(entry.getKey())) {
          header.mPathAttributes.put(entry.getKey(), entry.getValue());
        }
      }
    }
    final MutableGraph graph = factory.makeGraph(header.mContigOverlap, header.mContigAttributes, header.mPathAttributes);
    final Map<Long, Long> contigTranslationMap = readContigs(graph, md5sums);
    readPaths(graph, contigTranslationMap, md5sums);
    return graph;
  }

  /**
   * Read a graph into memory.
   * @param graphDir directory containing the graph
   * @return the graph
   * @throws NullPointerException if <code>graphDir</code> is null.
   * @throws IOException if an I/O error occurs.
   */
  public static Graph read(final StoreDirectory graphDir) throws IOException {
    return new GraphReader(graphDir).read(GraphFactory.DEFAULT, null, null);
  }

  /**
   * Read a graph into memory.
   * @param factory for constructing the new graph.
   * @param graphDir directory containing the graph
   * @return the graph
   * @throws NullPointerException if <code>graphDir</code> is null.
   * @throws IOException if an I/O error occurs.
   */
  public static Graph read(final GraphFactory factory, final StoreDirectory graphDir) throws IOException {
    return new GraphReader(graphDir).read(factory, null, null);
  }

  /**
   * Retrieve the guid from the graph
   * @param directory a directory containing the graph
   * @return a UUID read from the header file
   * @throws IOException when the UUID can't be read
   */
  public static UUID getUUID(StoreDirectory directory) throws IOException {
    final String headerFile;
    if (directory.childExists(HEADER_LEGACY)) {
      headerFile = HEADER_LEGACY;
    } else {
      headerFile = HEADER;
    }
    final StoreFile child = directory.child(headerFile);
    String uuidString = null;
    try (final BufferedReader reader = new BufferedReader(new InputStreamReader(child.inputStream()))) {
      String line;
      while ((line = reader.readLine()) != null) {
        final String[] split = StringUtils.split(line, '\t');
        if (split.length == 2 && "guid".equals(split[0])) {
          uuidString = split[1];
          break;
        }
      }
    }
    if (uuidString == null) {
      throw new IOException("UUID for graph is missing");
    }
    return UUID.fromString(uuidString);
  }

  /**
   * read a graph including additional attributes in the graph structure
   * @param factory for constructing the new graph.
   * @param graphDir directory containing the graph
   * @param additionalContigAttributes any contig attributes to add to those already in the graph
   * @param additionalPathAttributes any path attributes to add to those already in the graph
   * @return the graph
   * @throws NullPointerException if <code>graphDir</code> is null.
   * @throws IOException if an I/O error occurs.
   */
  public static Graph read(final GraphFactory factory
      , final StoreDirectory graphDir
      , Map<String, String> additionalContigAttributes
      , Map<String, String> additionalPathAttributes) throws IOException {
    return new GraphReader(graphDir).read(factory, additionalContigAttributes, additionalPathAttributes);
  }
}
