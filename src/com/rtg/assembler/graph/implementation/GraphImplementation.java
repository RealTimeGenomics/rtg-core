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

package com.rtg.assembler.graph.implementation;

import java.util.LinkedHashMap;
import java.util.Map;

import com.rtg.assembler.graph.Contig;
import com.rtg.assembler.graph.MutableGraph;
import com.rtg.assembler.graph.Path;
import com.rtg.assembler.graph.PathsIterator;
import com.rtg.util.Pair;
import com.rtg.util.array.ExtensibleIndex;
import com.rtg.util.array.bitindex.BitIndex;
import com.rtg.util.array.intindex.IntChunks;
import com.rtg.util.array.longindex.LongChunks;
import com.rtg.util.integrity.Exam;
import com.rtg.util.integrity.IntegralAbstract;
import com.rtg.variant.dna.DNARange;

/**
 * Uses index arrays. Not very object oriented but trades this with memory and speed efficiency.
 * The following diagram shows the relationship between the different tricky parts of the
 * data structures (doesn't include some parts like the deletions or the attributes).
 * <img src="doc-files/GraphImplementationDS/Slide1.jpg" />
 */
public class GraphImplementation extends IntegralAbstract implements MutableGraph {

  final int mContigOverlap;
  /**
   * For contig (a positive contig identifier) <code> mContigs[contig - 1]</code> is the first (0 based) position in
   * <code>mNt</code> containing a nucleotide for the contig. <code>mContigs[numberContigs()]</code> is defined and is one past the
   * last nucleotide stored in <code>mNt</code>.
   * As a consequence the 0 position is  0.
   */
  final ExtensibleIndex mContigs = new LongChunks(1);

  private final ExtensibleIndex mContigDeleted = new BitIndex(1, 1);

  /*
   * This is the same length as mContigs. The 0 position is not used and is 0.
   * For contig (a positive contig identifier) <code>mLastPath[contig]</code> points into <code>mPathIndex, mPathIndexOffset and mPathIndexPrev</code>.
   * All the paths for the contig can be found by following the chain of entries in <code>mPathIndexPrev</code> starting at <code>mLastPath[contig]</code>.
   * -1 is the null pointer (when a contig is not a member of any path).
   */
  final ExtensibleIndex mLastPath = new LongChunks(1);

  /** The concatenated nucleotides for all the contigs. See <code>mContigs</code>. */
  final ExtensibleIndex mNt = new BitIndex(0, 4);

  /**
   * For path (a positive path identifier) <code>mPath[path - 1]</code> is a (0 based) position in
   * <code>mPathIndex</code> (see there for details). <code>mPath[numberPaths()]</code> is defined and is one past the
   * last entry stored in <code>mPathIndex</code>.
   * As a consequence the 0 position is 0.
   */
  private final ExtensibleIndex mPath = new LongChunks(1);

  private final ExtensibleIndex mPathDeleted = new BitIndex(1, 1);


  /**
   * The position pointed at by <code>mPath[pathId - 1]</code> contains the <code>pathId</code> itself and is followed by one entry for each
   * contig in the path (there are always at least 2).
   */
  final ExtensibleIndex mPathIndex = new LongChunks(0);

  /**
   * This is in parallel with <code>mPathIndex</code>.
   * Given <code>x = mPath[pathId -1]</code>.
   * For 0 &le; i &le; path length, <code>mPathIndexOffset[x + i] == i</code>.
   */
  final ExtensibleIndex mPathIndexOffset = new IntChunks(0);

  /**
   * This is in parallel with <code>mPathIndex</code>.
   * Connects all the index entries for the paths that contain a particular contig.
   * -1 is the null pointer.
   * Given <code>x = mPath[pathId -1], mPathIndexPrev[x] == -1</code>.
   */
  final ExtensibleIndex mPathIndexPrev = new LongChunks(0);

  private final Map<String, String> mContigAttributes;

  private final Map<Pair<String, Long>, String> mContigAttributeMap = new LinkedHashMap<>();

  private final Map<String, String> mPathAttributes;

  private final Map<Pair<String, Long>, String> mPathAttributeMap = new LinkedHashMap<>();

  /**
   * @param contigOverlap the number of bases adjacent contigs overlap
   * @param contigAttributes key value pairs for contig attributes.
   * @param pathAttributes key value pairs for path attributes.
   */
  public GraphImplementation(int contigOverlap, final Map<String, String> contigAttributes, final Map<String, String> pathAttributes) {
    mContigOverlap = contigOverlap;
    mContigAttributes = new LinkedHashMap<>();
    mContigAttributes.putAll(contigAttributes);
    mPathAttributes = new LinkedHashMap<>();
    mPathAttributes.putAll(pathAttributes);
    mPathAttributes.put("deleted", "true if this path has been deleted");
    mContigAttributes.put("deleted", "true if this contig has been deleted");
    mLastPath.set(0, -1);
  }

  /**
   * Check that a signed contig is valid and return its absolute value.
   * @param contigId a putative contig identifier.
   * @return the positive contig identifier.
   */
  final long absContig(long contigId) {
    final long abs;
    if (contigId < 0) {
      abs = -contigId;
    } else if (contigId > 0) {
      abs = contigId;
    } else {
      throw new IllegalArgumentException();
    }
    if (abs >= mContigs.length()) {
      throw new IllegalArgumentException(abs + ":" + mContigs.length());
    }
    return abs;
  }

  /**
   * Check that a signed contig is valid and return its absolute value.
   * @param pathId a putative contig identifier.
   * @return the positive contig identifier.
   */
  final long absPath(long pathId) {
    final long abs;
    if (pathId < 0) {
      abs = -pathId;
    } else if (pathId > 0) {
      abs = pathId;
    } else {
      throw new IllegalArgumentException();
    }
    if (abs >= mPath.length()) {
      throw new IllegalArgumentException(abs + ":" + mPath.length());
    }
    return abs;
  }

  @Override
  public int contigOverlap() {
    return mContigOverlap;
  }

  @Override
  public final long numberContigs() {
    return mContigs.length() - 1;
  }

  @Override
  public final long numberPaths() {
    return mPath.length() - 1;
  }

  @Override
  public final boolean contigDeleted(long contigId) {
    final long acontig = absContig(contigId);
    final long del = mContigDeleted.get(acontig);
    assert del == 0 || del == 1;
    return del == 1;
  }

  @Override
  public final int contigLength(long contigId) {
    final long con = absContig(contigId);
    final long end = mContigs.get(con);
    final long start = mContigs.get(con - 1);
    final long length = end - start;
    //assert length <= Integer.MAX_VALUE;
    return (int) length;
  }

  @Override
  public final byte nt(long contigId, int index) {
    if (index < 0) {
      throw new IllegalArgumentException(contigId + ":" + index);
    }
    final byte nt;
    if (contigId > 0) {
      final long start = mContigs.get(contigId - 1);
      final long end = mContigs.get(contigId);
      final long ix = start + index;
      if (ix >= end) {
        throw new IllegalArgumentException(contigId + ":" + index);
      }
      nt = (byte) mNt.get(ix);
    } else if (contigId < 0) {
      final long end = mContigs.get(-contigId) - 1;
      final long start = mContigs.get(-contigId - 1) - 1;
      final long ix = end - index;
      if (ix <= start) {
        throw new IllegalArgumentException(contigId + ":" + index);
      }
      final byte ntc = (byte) mNt.get(ix);
      nt = DNARange.complement(ntc);
    } else {
      throw new IllegalArgumentException();
    }
    return nt;
  }

  @Override
  public final Contig contig(long contigId) {
    if (contigId < 0) {
      return new ContigReverse(this, contigId);
    } else if (contigId > 0) {
      return new ContigForward(this, contigId);
    } else {
      throw new IllegalArgumentException();
    }
  }

  @Override
  public final Map<String, String> contigAttributes() {
    return mContigAttributes;
  }

  @Override
  public String contigAttribute(long contigId, String attribute) {
    if (!mContigAttributes.containsKey(attribute)) {
      throw new IllegalArgumentException();
    }
    if (attribute.equals("deleted")) {
      return contigDeleted(contigId) ? Boolean.toString(true) : null;
    }
    final long acontig = absContig(contigId);
    return mContigAttributeMap.get(new Pair<>(attribute, acontig));
  }

  @Override
  public final Map<String, String> pathAttributes() {
    return mPathAttributes;
  }

  @Override
  public String pathAttribute(long pathId, String attribute) {
    if (!mPathAttributes.containsKey(attribute)) {
      throw new IllegalArgumentException();
    }
    if (attribute.equals("deleted")) {
      return pathDeleted(pathId) ? Boolean.toString(true) : null;
    }
    final long apath = absPath(pathId);
    return mPathAttributeMap.get(new Pair<>(attribute, apath));
  }

  @Override
  public final boolean pathDeleted(long pathId) {
    final long apath = absPath(pathId);
    final long del = mPathDeleted.get(apath);
    final boolean delb = del == 1;
    assert del == 0 || del == 1;
    return delb;
  }

  @Override
  public final PathsIterator iterator() {
    return new IteratorLocal(this);
  }

  @Override
  public final PathsIterator paths(long contigId) {
    final IteratorLocal it = new IteratorLocal(this);
    it.set(contigId);
    return it;
  }

  @Override
  public final PathsIterator paths(long contigId, boolean showDeleted) {
    final IteratorLocal it = new IteratorLocal(this);
    it.set(contigId, showDeleted);
    return it;
  }

  @Override
  public final int pathLength(long pathId) {
    final long apath = absPath(pathId);
    final long length = mPath.get(apath) - mPath.get(apath - 1) - 1;
    //assert length <= Integer.MAX_VALUE;
    return (int) length;
  }

  @Override
  public final long pathContig(long pathId, int index) {
    if (index < 0) {
      throw new IllegalArgumentException();
    }
    final long contig;
    if (pathId > 0) {
      final long start = mPath.get(pathId - 1);
      final long end = mPath.get(pathId);

      final long ix = start + 1 + index;
      if (ix >= end) {
        throw new IllegalArgumentException(pathId + ":" + index);
      }
      contig = mPathIndex.get(ix);
    } else if (pathId < 0) {
      final long end = mPath.get(-pathId) - 1;
      final long start = mPath.get(-pathId - 1);
      final long ix = end - index;
      if (ix <= start) {
        throw new IllegalArgumentException(pathId + ":" + index);
      }
      contig = -mPathIndex.get(ix);
      //System.err.println("end=" + end + " contig=" + contig);
    } else {
      throw new IllegalArgumentException();
    }
    //System.err.println("pathId=" + pathId + " index=" + index + " contig=" + contig);
    return contig;
  }

  @Override
  public final Path path(long pathId) {
    return new PathLocal(this, pathId);
  }

  //The methods from MutableGraph follow

  @Override
  public final long addContig(Contig contig) {
    final long next = mContigs.length();
    final long start = mContigs.get(next - 1);
    final long end = start + contig.length();
    contigExpand();
    mNt.extendBy(contig.length());
    for (int i = 0; i < contig.length(); i++) {
      mNt.set(start + i, contig.nt(i));
    }
    mContigs.set(next, end);
    mLastPath.set(next, -1);
    assert mNt.length() == end;
    return next;
  }

  /**
   * Expand contig and any indexes that parallel it by 1 entry.
   */
  protected void contigExpand() {
    mContigs.extendBy(1);
    mLastPath.extendBy(1);
    mContigDeleted.extendBy(1);
  }

  @Override
  public final void deleteContig(long contigId) throws IllegalArgumentException, IllegalStateException {
    final long acontig = absContig(contigId);
    mContigDeleted.set(acontig, 1);
    final PathsIterator it = paths(contigId, false);
    while (true) {
      final long path = it.nextPathId();
      if (path == 0) {
        break;
      }
      deletePath(path);
    }
  }

  @Override
  public void setContigAttribute(long contigId, String attribute, String value) throws IllegalArgumentException, IllegalStateException {
    if (!mContigAttributes.containsKey(attribute)) {
      throw new IllegalArgumentException();
    }
    if (attribute.equals("deleted")) {
      if (Boolean.parseBoolean(value)) {
        deleteContig(contigId);
      }
    }

    final long acontig = absContig(contigId);
    mContigAttributeMap.put(new Pair<>(attribute, acontig), value);
  }

  @Override
  public final long addPath(Path path) {
    final long next = mPath.length();
    final long start = mPath.get(next - 1);
    final int incr = path.length() + 1;
    final long end = start + incr;
    pathExpand();
    pathIndexExpand(incr);
    mPath.set(next, end);
    mPathIndex.set(start, next);
    mPathIndexOffset.set(start, 0);
    mPathIndexPrev.set(start, -1);
    for (int i = 0; i < path.length(); i++) {
      final long index = start + i + 1;
      final long contig = path.contig(i);
      final long acontig = absContig(contig);
      mPathIndex.set(index, contig);
      mPathIndexOffset.set(index, i + 1);
      final long prev = mLastPath.get(acontig);
      mPathIndexPrev.set(index, prev);
      mLastPath.set(acontig, index);
    }
    return next;
  }

  protected void pathExpand() {
    mPath.extendBy(1);
    mPathDeleted.extendBy(1);
  }

  protected void pathIndexExpand(final long incr) {
    mPathIndex.extendBy(incr);
    mPathIndexOffset.extendBy(incr);
    mPathIndexPrev.extendBy(incr);
  }

  @Override
  public final void deletePath(long pathId) throws IllegalArgumentException, IllegalStateException {
    final long apath = absPath(pathId);
    mPathDeleted.set(apath, 1);
  }

  @Override
  public void setPathAttribute(long pathId, String attribute, String value) {
    if (!mPathAttributes.containsKey(attribute)) {
      throw new IllegalArgumentException();
    }
    if (attribute.equals("deleted")) {
      if (Boolean.parseBoolean(value)) {
        deletePath(pathId);
      }
    }
    final long acontig = absPath(pathId);
    mPathAttributeMap.put(new Pair<>(attribute, acontig), value);
  }

  @Override
  public void addPathAttribute(String name, String description) {
    mPathAttributes.put(name, description);
  }

  @Override
  public void addContigAttribute(String name, String description) {
    mContigAttributes.put(name, description);
  }

  @Override
  public final MutableGraph compact() {
    final MutableGraph newGraph = new GraphImplementation(mContigOverlap, contigAttributes(), pathAttributes());
    compact(newGraph);
    return newGraph;
  }

  @Override
  public final void compact(MutableGraph newGraph) {
    assert newGraph.numberContigs() == 0;
    final ExtensibleIndex transx = new LongChunks(numberContigs() + 1);
    long count = 1;
    for (long l = 1; l <= numberContigs(); l++) {
      if (!contigDeleted(l)) {
        newGraph.addContig(contig(l));
        transx.set(l, count);
        for (final String attribute : contigAttributes().keySet()) {
          final String value = contigAttribute(l, attribute);
          if (value != null) {
            newGraph.setContigAttribute(count, attribute, value);
          }
        }
        count++;
      }
    }
    for (long l = 1; l <= numberPaths(); l++) {
      if (!pathDeleted(l)) {
        final Path path = path(l);
        final long[] contigs = new long[path.length()];
        for (int i = 0; i < contigs.length; i++) {
          final long contig = path.contig(i);
          assert contig != 0;
          final long newContig;
          if (contig > 0) {
            newContig = transx.get(contig);
          } else {
            newContig = -transx.get(-contig);
          }
          contigs[i] = newContig;
        }
        final Path newPath = new PathArray(contigs);
        final long newPathId = newGraph.addPath(newPath);
        for (final String attribute : pathAttributes().keySet()) {
          final String value = pathAttribute(l, attribute);
          if (value != null) {
            newGraph.setPathAttribute(newPathId, attribute, value);
          }
        }
      }
    }
  }

  @Override
  public boolean globalIntegrity() {
    integrity();

    //check the indexes independently to make sure they all have reasonable values
    //for each entry in them.
    long last0 = 0;
    final long contigsLength = mContigs.length();
    for (long l = 1; l < contigsLength; l++) {
      final long next = mContigs.get(l);
      Exam.assertTrue(next >= last0);
      last0 = next;
    }

    for (long l = 0; l < mNt.length(); l++) {
      final long nt = mNt.get(l);
      DNARange.RANGE.check(nt);
    }

    long last1 = 0;
    final long pathLength = mPath.length();
    for (long l = 1; l < pathLength; l++) {
      final long next = mPath.get(l);
      Exam.assertTrue(next >= last1);
      last1 = next;
    }

    final long pathIndexLength = mPathIndex.length();
    for (long l = 0; l < mLastPath.length(); l++) {
      final long index = mLastPath.get(l);
      Exam.assertTrue(l + ":" + index, index >= -1 && index < pathIndexLength);
    }

    long lastOffset = 2;
    for (long l = 0; l < pathIndexLength; l++) {
      final long offset = mPathIndexOffset.get(l);
      final long index = mPathIndex.get(l);

      Exam.assertTrue(offset >= 0);
      if (offset == 0) {
        Exam.assertTrue(1 <= index && index < pathLength);
        Exam.assertTrue(lastOffset >= 2);
      } else {
        final long aindex = absContig(index);
        Exam.assertTrue(aindex + ":" + contigsLength, 1 <= aindex && aindex < contigsLength);
        Exam.assertEquals(lastOffset, offset - 1);
      }
      lastOffset = offset;
    }

    for (long l = 0; l < pathIndexLength; l++) {
      final long prev = mPathIndexPrev.get(l);
      if (prev != -1) {
        Exam.assertTrue(0 <= prev && prev < l);
      }
    }

    //check linkages between indexes
    //lastPath and pathIndex have bidirectional pointers.
    for (long l = 1; l < mLastPath.length(); l++) {
      final long index = mLastPath.get(l);
      if (index == -1) {
        continue;
      }
      final long back = mPathIndex.get(index);
      final long aback = absContig(back);
      Exam.assertEquals(l, aback);
    }

    //path and pathIndex have bidirectional pointers.
    for (long l = 0; l < pathLength - 1; l++) {
      final long path = mPath.get(l);
      final long back = absContig(mPathIndex.get(path));
      Exam.assertEquals("" + l, l + 1, back);
      final long offset = mPathIndexOffset.get(path);
      Exam.assertEquals(0, offset);
    }

    //The chains of pathIndexPrev all point back at the same contig
    for (long l = 1; l < mLastPath.length(); l++) {
      long ix = mLastPath.get(l);
      while (ix != -1) {
        final long back = mPathIndex.get(ix);
        final long aback = absContig(back);
        Exam.assertEquals(l, aback);
        ix = mPathIndexPrev.get(ix);
      }
    }

    return true;
  }

  /**
   * @return the sum of bytes used in indexes
   */
  public long bytes() {
    return mContigs.bytes() + mPath.bytes() + mContigDeleted.bytes() + mPathDeleted.bytes() + mPathIndex.bytes() + mNt.bytes() + mLastPath.bytes() + mPathIndexOffset.bytes() + mPathIndexPrev.bytes();
  }

  @Override
  public boolean integrity() {
    Exam.assertTrue(mContigs.length() >= 1);
    Exam.assertEquals(0, mContigs.get(0));
    Exam.assertEquals(mContigs.get(mContigs.length() - 1), mNt.length());
    Exam.assertEquals(mContigs.length(), mLastPath.length());
    Exam.assertEquals(-1, mLastPath.get(0));
    Exam.assertEquals(mContigDeleted.length(), mLastPath.length());

    Exam.assertTrue(mPath.length() >= 1);
    Exam.assertEquals(mPathDeleted.length(), mPath.length());
    Exam.assertEquals(0, mPath.get(0));
    Exam.assertEquals(mPath.get(mPath.length() - 1), mPathIndex.length());

    Exam.assertTrue(mPathIndex.length() >= 2 * (mPath.length() - 1));
    Exam.assertEquals(mPathIndex.length(), mPathIndexOffset.length());
    if (mPathIndexOffset.length() > 0) {
      Exam.assertEquals(0, mPathIndexOffset.get(0));
    }
    Exam.assertEquals(mPathIndex.length(), mPathIndexPrev.length());
    return true;
  }
}
