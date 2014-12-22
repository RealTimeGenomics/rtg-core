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

package com.rtg.variant.eval;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Stack;
import java.util.TreeMap;
import java.util.TreeSet;

import com.rtg.util.BasicLinkedListNode;
import com.rtg.util.Pair;
import com.rtg.util.Utils;
import com.rtg.util.diagnostic.Diagnostic;
import com.rtg.util.diagnostic.NoTalkbackSlimException;
import com.rtg.util.integrity.IntegralAbstract;

/**
 * Record a path that reconciles two sequences of variants. Also provides
 * method to find a best path.
 *
 */
public class Path extends IntegralAbstract implements Comparable<Path> {
  private static final int MAX_COMPLEXITY = 5000;

  private final HalfPath mCalledPath;
  private final HalfPath mBaselinePath;
  //private final byte[] mTemplate;

  private final BasicLinkedListNode<Integer> mSyncPointList;

  static class SyncPoint implements Comparable<SyncPoint> {
    private final int mPos;
    private final double mCalledTPCount;
    private final double mBaselineTPCount;

    public SyncPoint(int pos, int calledCounts, int baselineCounts) {
      mPos = pos;
      mCalledTPCount = calledCounts;
      mBaselineTPCount = baselineCounts;
    }

    int getPos() {
      return mPos;
    }

    @Override
    public String toString() {
      return "SyncPoint [mPos=" + mPos + ", mCalledTPCount=" + mCalledTPCount + ", mBaselineTPCount=" + mBaselineTPCount + "]";
    }

    @Override
    public int hashCode() {
      return mPos;
    }

    @Override
    public boolean equals(Object obj) {
      if (obj == null) {
        return false;
      }
      if (!(obj instanceof SyncPoint)) {
        return false;
      }
      return this.compareTo((SyncPoint) obj) == 0;
    }

    @Override
    public int compareTo(SyncPoint o) {
      if (this.mPos < o.mPos) {
        return -1;
      } else if (this.mPos > o.mPos) {
        return 1;
      }
      return 0;
    }

  }

  /**
   * Construct a path object for the specified template
   *
   * @param template the template sequence to build a path for.
   */
  public Path(byte[] template) {
    //mTemplate = template;
    mCalledPath = new HalfPath(template);
    mBaselinePath = new HalfPath(template);
    mSyncPointList = null;
  }

  private Path(Path parent, BasicLinkedListNode<Integer> syncPoints) {
    mCalledPath = new HalfPath(parent.mCalledPath);
    mBaselinePath = new HalfPath(parent.mBaselinePath);
    //mTemplate = parent.mTemplate;
    mSyncPointList = syncPoints;
  }

  /**
   * @return true if the path has finished
   */
  public boolean finished() {
    return mCalledPath.finished() && mBaselinePath.finished();
  }

  /**
   * @return true if all positions are the same
   */
  public boolean inSync() {
    if (mCalledPath.compareHaplotypePositions() != 0) {
      return false;
    }
    if (mBaselinePath.compareHaplotypePositions() != 0) {
      return false;
    }
    if (mCalledPath.getPosition() != mBaselinePath.getPosition()) {
      return false;
    }
    if (mCalledPath.getPosition() <= mCalledPath.getVariantEndPosition()) {
      return false;
    }

    if (mBaselinePath.getPosition() <= mBaselinePath.getVariantEndPosition()) {
      return false;
    }
    return true;
  }

  static Path better(Path first, Path second) {
    final BasicLinkedListNode<OrientedVariant> firstIncluded =  first == null ? null : first.mCalledPath.getIncluded();
    final BasicLinkedListNode<OrientedVariant> secondIncluded = second == null ? null : second.mCalledPath.getIncluded();
    final int firstSize = firstIncluded == null ? 0 : firstIncluded.size();
    final int secondSize = secondIncluded == null ? 0 : secondIncluded.size();
    if (firstSize == secondSize) {
      final BasicLinkedListNode<OrientedVariant> firstBaseline =  first == null ? null : first.mBaselinePath.getIncluded();
      final BasicLinkedListNode<OrientedVariant> secondBaseline = second == null ? null : second.mBaselinePath.getIncluded();
      final int firstBaseSize = firstBaseline == null ? 0 : firstBaseline.size();
      final int secondBaseSize = secondBaseline == null ? 0 : secondBaseline.size();
      if (firstBaseSize == secondBaseSize) {
        if (firstBaseline != null && secondBaseline != null) {
          if (firstBaseline.getValue().isAlleleA() && !secondBaseline.getValue().isAlleleA()) {
            return first;
          } else if (secondBaseline.getValue().isAlleleA()) {
            return second;
          }
        }
      }
      return firstBaseSize > secondBaseSize ? first : second;
    }
    return firstSize > secondSize ? first : second;
  }

  /**
   * Find the path through the two sequences of variants that best reconciles
   * them.
   *
   * @param template original reference sequence.
   * @param templateName name of the current reference sequence
   * @param calledVariants first sequence of variants to be applied to the template.
   * @param baseLineVariants second sequence of variants to be applied to the template.
   * @param <T> the type parameter
   * @return the best path (non-null).
   */
  public static <T extends Variant> Path bestPath(byte[] template, String templateName, Collection<T> calledVariants, Collection<T> baseLineVariants) {
    // make it easy to find variations
    final TreeMap<Integer, Variant> calledMap = buildMap(calledVariants);
    final TreeMap<Integer, Variant> baselineMap = buildMap(baseLineVariants);
    final TreeSet<Path> sortedPaths = new TreeSet<>();
    sortedPaths.add(new Path(template));
    Path best = null;
    int maxPaths = 0;
    String maxPathsRegion = "";
    int currentMax = 0;
    int lastSyncPos = 0;
    while (sortedPaths.size() > 0) {
      currentMax = Math.max(currentMax, sortedPaths.size());
      final Path head = sortedPaths.pollFirst();
      if (sortedPaths.size() == 0) { // Only one path currently in play
        final int currentSyncPos = head.mCalledPath.getPosition();
        if (currentMax > maxPaths) {
          maxPathsRegion = templateName + ":" + (lastSyncPos + 1) + "-" + (currentSyncPos + 1);
          maxPaths = currentMax;
          Diagnostic.developerLog("Maximum path complexity now " + maxPaths + ", at " + maxPathsRegion);
        }
        currentMax = 0;
        lastSyncPos = currentSyncPos;
      } else if (sortedPaths.size() > MAX_COMPLEXITY) {
        final int currentMaxPos = Math.max(head.mCalledPath.getPosition(), head.mBaselinePath.getPosition());
        throw new NoTalkbackSlimException("Evaluation got far too hard around reference sequence " + templateName + ":" + (lastSyncPos + 1) + "-" + (currentMaxPos + 1) + ". Look at the variants in (and after) this region, and try filtering these out.");
      }
      if (head.finished()) {
        // Path is done. Remember the best
        final BasicLinkedListNode<Integer> syncPoints = new BasicLinkedListNode<>(head.mCalledPath.getPosition(), head.mSyncPointList);
        best = better(best, new Path(head, syncPoints));
        continue;
      }
      final Variant aVar = nextVariant(head.mCalledPath, calledMap);
      if (aVar != null && head.mCalledPath.isNew(aVar)) {
        //Adding a new variant to A side
        addIfBetter(head.addAVariant(aVar), sortedPaths);
        continue;
      }
      final Variant bVar = nextVariant(head.mBaselinePath, baselineMap);
      if (bVar != null && head.mBaselinePath.isNew(bVar)) {
        //Adding a new variant to B side
        addIfBetter(head.addBVariant(bVar), sortedPaths);
        continue;
      }

      head.step();

      if (head.inSync()) {
        skipToNextVariant(template, calledMap, baselineMap, head);
      }
      
      if (head.matches()) {
        addIfBetter(head, sortedPaths);
      }
    }
    Diagnostic.userLog("Reference " + templateName + " had maximum path complexity of " + maxPaths + " at " + maxPathsRegion);
    return best;
  }

  /**
   * Move the path to just before the next variant
   * @param template the template to skip along
   * @param aMap list of variants for the A haplotype
   * @param bMap list of variants for the B haplotype
   * @param head the path to skip forward
   */
  private static void skipToNextVariant(byte[] template, final TreeMap<Integer, Variant> aMap, final TreeMap<Integer, Variant> bMap, final Path head) {
    final Variant aNext = futureVariant(head.mCalledPath, aMap);
    final Variant bNext = futureVariant(head.mBaselinePath, bMap);
    final int lastTemplatePos = template.length - 1;
    // -1 because we want to be before the position
    final int nextPos = Math.min(Math.min(aNext != null ? aNext.getStart() : lastTemplatePos, bNext != null ? bNext.getStart() : lastTemplatePos), lastTemplatePos) - 1;
    if (nextPos > head.mCalledPath.getPosition()) {
      head.moveForward(nextPos);
    }
  }

  static Variant nextVariant(HalfPath path, Map<Integer, Variant> map) {
    return map.get(Math.max(path.getVariantEndPosition(), path.getPosition() + 1));
  }
  static Variant futureVariant(HalfPath path, TreeMap<Integer, Variant> map) {
    final Entry<Integer, Variant> entry = map.ceilingEntry(path.getPosition());
    return entry == null ? null : entry.getValue();
  }

  private static void addIfBetter(Collection<Path> add, TreeSet<Path> sortedPaths) {
        for (final Path p : add) {
          addIfBetter(p, sortedPaths);
        }
  }

  private static void addIfBetter(Path add, TreeSet<Path> sortedPaths) {
    if (sortedPaths.contains(add)) {
      final Path other = sortedPaths.floor(add);
      sortedPaths.remove(add);
      sortedPaths.add(better(add, other));
    } else {
      sortedPaths.add(add);
    }
  }

  static <T extends Variant> TreeMap<Integer, Variant> buildMap(Collection<T> variations) {
    final TreeMap<Integer, Variant> map = new TreeMap<>();
    for (final Variant v : variations) {
      // TODO when you have an insert immediately prior to a snp/mnp, they end up with the same start position,
      // but don't trigger the overlapping code during variant loading.
      // This means that you don't get a warning but the snp/mnp is the only variant that ends up in the map
      // due to the map.put replacing the existing entry. For now log these to see how often they occur.
      final Variant oldV = map.put(v.getStart(), v);
      if (oldV != null) {
        Diagnostic.developerLog("Variant was bumped due to another variant starting at the same position.\nCurrent variant: " + v + "\nBumped variant:  " + oldV);
      }
    }
    return map;
  }

  @Override
  public void toString(StringBuilder sb) {
  }

  @Override
  public boolean integrity() {
    return true;
  }
  private static List<OrientedVariant> getIncluded(HalfPath path) {
    final List<OrientedVariant> list = new ArrayList<>();
    if (path.getIncluded() == null) {
      return list;
    }
    for (final OrientedVariant orientedVariant : path.getIncluded()) {
      list.add(0, orientedVariant);
    }
    return list;
  }

  private static List<Variant> getExcluded(HalfPath path) {
    final List<Variant> list = new ArrayList<>();
    if (path.getExcluded() == null) {
      return list;
    }
    for (final Variant variant : path.getExcluded()) {
      list.add(0, variant);
    }
    return list;
  }

  /**
   * @return the set of variants included on A side
   */
  public List<OrientedVariant> getCalledIncluded() {
    return getIncluded(mCalledPath);
  }

  /**
   * @return the set of variants excluded on A side
   */
  public List<Variant> getCalledExcluded() {
    return getExcluded(mCalledPath);
  }

  /**
   * @return the set of variants included on B side
   */
  public List<OrientedVariant> getBaselineIncluded() {
    return getIncluded(mBaselinePath);
  }

  /**
   * @return the set of variants excluded on B side
   */
  public List<Variant> getBaselineExcluded() {
    return getExcluded(mBaselinePath);
  }

  void include(boolean side, OrientedVariant var) {
    if (side) {
      mCalledPath.include(var);
    } else {
      mBaselinePath.include(var);
    }
  }

  void exclude(boolean side, Variant var) {
    if (side) {
      mCalledPath.exclude(var);
    } else {
      mBaselinePath.exclude(var);
    }
  }

  Collection<Path> addVariant(boolean side, Variant var) {
    final ArrayList<Path> paths = new ArrayList<>();
    final BasicLinkedListNode<Integer> syncPoints;
    if (this.inSync()) {
      syncPoints = new BasicLinkedListNode<>(this.mCalledPath.getPosition(), this.mSyncPointList);
    } else {
      syncPoints = this.mSyncPointList;
    }
    final Path exclude = new Path(this, syncPoints);
    exclude.exclude(side, var);
    paths.add(exclude);
    final Path include = new Path(this, syncPoints);
    include.include(side, new OrientedVariant(var, true));
    assert !include.equals(exclude);
    assert include.compareTo(exclude) != 0;
    paths.add(include);
    if (var.ntAlleleB() != null) {
      final Path hetero = new Path(this, syncPoints);
      hetero.include(side, new OrientedVariant(var, false));
      paths.add(hetero);
    }
    return paths;
  }

  Collection<Path> addAVariant(Variant var) {
    return addVariant(true, var);
  }

  Collection<Path> addBVariant(Variant var) {
    return addVariant(false, var);
  }
  void moveForward(int position) {
    mCalledPath.moveForward(position);
    mBaselinePath.moveForward(position);
  }

  void step() {
    if (mCalledPath.compareHaplotypePositions() > 0) {
      // make B haplotype catch up to A
      minusStep();
    } else if (mCalledPath.compareHaplotypePositions() < 0) {
      // make A haplotype catch up to B
      plusStep();
    } else {
      // step both
      mCalledPath.step();
      mBaselinePath.step();
    }
  }

  void minusStep() {
    mCalledPath.haplotypeBStep();
    mBaselinePath.haplotypeBStep();
  }

  void plusStep() {
    mCalledPath.haplotypeAStep();
    mBaselinePath.haplotypeAStep();
  }

  boolean matches() {
    if (!(mCalledPath.finishedHaplotypeA() || mBaselinePath.finishedHaplotypeA()) && mCalledPath.nextHaplotypeABase() != mBaselinePath.nextHaplotypeABase()) {
      return false;
    }
    if (!(mCalledPath.finishedHaplotypeB() || mBaselinePath.finishedHaplotypeB()) && mCalledPath.nextHaplotypeBBase() != mBaselinePath.nextHaplotypeBBase()) {
      return false;
    }
    return true;
  }

  @Override
  public int compareTo(Path o) {
    final int aComp = mCalledPath.compareTo(o.mCalledPath);
    if (aComp != 0) {
      return aComp;
    }
    return mBaselinePath.compareTo(o.mBaselinePath);
  }

  @Override
  public boolean equals(Object o) {
    if (o == null) {
      return false;
    }

    if (!o.getClass().equals(getClass())) {
      return false;
    }
    return compareTo((Path) o) == 0;
  }

  @Override
  public int hashCode() {
    return Utils.hash(new Object[] {mCalledPath, mBaselinePath});
  }

  @Override
  public String toString() {

    return "Path: called=" + mCalledPath + " baseline=" + mBaselinePath + HalfPath.varListToString(this.mSyncPointList);
  }
  public List<SyncPoint> getSyncPoints() {
    return getSyncPointsList(mSyncPointList, getCalledIncluded(), getBaselineIncluded());
  }

  // Counts the baseline and called variants between each sync point
  static List<SyncPoint> getSyncPointsList(final BasicLinkedListNode<Integer> syncpoints, final List<OrientedVariant> baseLine, final List<OrientedVariant> called) {
    // Push all the sync points onto a stack so we can process in reverse order.
    final Stack<Integer> stack = new Stack<>();
    for (final Integer syncpoint : syncpoints) {
      stack.push(syncpoint);
    }
    final List<SyncPoint> list = new ArrayList<>();
    int basePos = 0;
    int callPos = 0;
    while (!stack.isEmpty()) {
      final int loc = stack.pop();
      int baseLineCount = 0;
      int calledCount = 0;
      while (baseLine.size() > basePos && baseLine.get(basePos).variant().getStart() <= loc) {
        baseLineCount++;
        basePos++;
      }

      while (called.size() > callPos && called.get(callPos).variant().getStart() <= loc) {
        calledCount++;
        callPos++;
      }
      list.add(new SyncPoint(loc, calledCount, baseLineCount));
    }
    return list;
  }

  /**
   * Find a weighting for all the calls in a path.
   * this is done by sync points, within each <code>SyncPoint</code>
   *
   *       weight = number of TP in baseline / number of TP in called
   *
   *  this will assure that the total number of TP we output will always reflect number of TP in baseline file
   *
   *  return 2 lists first for TP and second for FP, we create FP calls when current call is included but does not correspond to a TP
   *
   * @param best best path
   * @param calledTruePositives called true positives
   * @param baseLineTruePositives baseline true positives
   *
   * @return pair of true positive and false positives
   */
  static Pair<List<OrientedVariant>, List<OrientedVariant>> calculateWeights(final Path best, final List<OrientedVariant> calledTruePositives, final List<OrientedVariant> baseLineTruePositives) {
    assert best.mSyncPointList.size() >= 1;
    final List<SyncPoint> syncpoints = getSyncPointsList(best.mSyncPointList, baseLineTruePositives, calledTruePositives);
    assert syncpoints.size() == best.mSyncPointList.size();

    final List<OrientedVariant> tp = new ArrayList<>();
    final List<OrientedVariant> fp = new ArrayList<>();

    final Iterator<SyncPoint> syncIterator = syncpoints.iterator();
    int syncStart = 0;
    SyncPoint syncpoint = syncIterator.next();
    for (final OrientedVariant v : calledTruePositives) {
      while (syncpoint.mPos < v.getStart() && syncIterator.hasNext()) { // Jump to sync point entry containing this variant
        syncStart = syncpoint.mPos;
        syncpoint = syncIterator.next();
      }
      if (syncpoint.mBaselineTPCount == 0) {
        Diagnostic.developerLog("Best path called variant was bumped to FP due to no baseline TP within sync region "
            + v.getSequenceName() + ":" + (syncStart + 1) + "-" + (syncpoint.mPos + 1) + "\n"
            + "Bumped variant:  " + v);
        fp.add(v);
      } else {
        v.setWeight(syncpoint.mBaselineTPCount / syncpoint.mCalledTPCount);
        tp.add(v);
      }

    }
    return new Pair<>(tp, fp);
  }




}
