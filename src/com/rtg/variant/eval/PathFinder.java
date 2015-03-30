/*
 * Copyright (c) 2015. Real Time Genomics Limited.
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

import java.util.Collection;
import java.util.Map;
import java.util.TreeMap;
import java.util.TreeSet;

import com.rtg.launcher.GlobalFlags;
import com.rtg.util.BasicLinkedListNode;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Find the path through the two sequences of variants that best reconciles them.
 */
public class PathFinder {

  private static final int MAX_COMPLEXITY = GlobalFlags.getIntegerValue(GlobalFlags.VCFEVAL_MAX_PATHS);

  private final byte[] mTemplate;
  private final String mTemplateName;
  private final TreeMap<Integer, Variant> mCalledVariants;
  private final TreeMap<Integer, Variant> mBaseLineVariants;

  <T extends Variant> PathFinder(byte[] template, String templateName, Collection<T> calledVariants, Collection<T> baseLineVariants) {
    mTemplate = template;
    mTemplateName = templateName;
    mCalledVariants = buildMap(calledVariants);
    mBaseLineVariants = buildMap(baseLineVariants);
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
    return new PathFinder(template, templateName, calledVariants, baseLineVariants).bestPath();
  }

  private Path bestPath() {
    // make it easy to find variations
    final TreeSet<Path> sortedPaths = new TreeSet<>();
    sortedPaths.add(new Path(mTemplate));
    Path best = null;
    int maxPaths = 0;
    String maxPathsRegion = "";
    int currentMax = 0;
    int currentMaxPos = 0;
    Path lastSyncPath = null;
    int lastSyncPos = 0;
    String lastWarnMessage = null;
    while (sortedPaths.size() > 0) {
      currentMax = Math.max(currentMax, sortedPaths.size());
      Path head = sortedPaths.pollFirst();
      //System.err.println("Size: " + (sortedPaths.size() + 1) + " Range:" + (lastSyncPos + 1) + "-" + (currentMaxPos + 1) + "\n\nHead: " + head);
      if (sortedPaths.size() == 0) { // Only one path currently in play
        if (lastWarnMessage != null) { // Issue a warning if we encountered problems during the previous region
          Diagnostic.warning(lastWarnMessage);
          lastWarnMessage = null;
        }
        final int currentSyncPos = head.mCalledPath.getPosition();
        if (currentMax > 1) {
          Diagnostic.developerLog("Most recent +5 region " + mTemplateName + ":" + (lastSyncPos + 1) + "-" + (currentSyncPos + 1));
        }
        if (currentMax > maxPaths) {
          maxPathsRegion = mTemplateName + ":" + (lastSyncPos + 1) + "-" + (currentSyncPos + 1);
          maxPaths = currentMax;
          Diagnostic.developerLog("Maximum path complexity now " + maxPaths + ", at " + maxPathsRegion);
        }
        currentMax = 0;
        lastSyncPos = currentSyncPos;
        lastSyncPath = head;
      } else if (sortedPaths.size() > MAX_COMPLEXITY) {
        lastWarnMessage = "Evaluation too complex (" + sortedPaths.size() + " unresolved paths) at reference region " + mTemplateName + ":" + (lastSyncPos + 1) + "-" + (currentMaxPos + 2) + ". Variants in this region will not be included in results.";
        sortedPaths.clear();    // Drop all paths currently in play
        head = lastSyncPath;    // Create new head containing path up until last sync point
        head.moveForward(currentMaxPos + 1);  // Skip to currentMaxPos
      }
      if (head.finished()) {
        // Path is done. Remember the best
        final BasicLinkedListNode<Integer> syncPoints = new BasicLinkedListNode<>(head.mCalledPath.getPosition(), head.mSyncPointList);
        best = better(best, new Path(head, syncPoints));
        continue;
      }
      final Variant aVar = nextVariant(head.mCalledPath, mCalledVariants);
      if (aVar != null && head.mCalledPath.isNew(aVar)) {
        currentMaxPos = Math.max(currentMaxPos, aVar.getStart());
        //Adding a new variant to A side
        //System.err.println("Add to A " + aVar);
        addIfBetter(head.addAVariant(aVar), sortedPaths);
        continue;
      }
      final Variant bVar = nextVariant(head.mBaselinePath, mBaseLineVariants);
      if (bVar != null && head.mBaselinePath.isNew(bVar)) {
        currentMaxPos = Math.max(currentMaxPos, bVar.getStart());
        //Adding a new variant to B side
        //System.err.println("Add to B " + bVar);
        addIfBetter(head.addBVariant(bVar), sortedPaths);
        continue;
      }

      head.step();

      if (head.inSync()) {
        skipToNextVariant(mTemplate, mCalledVariants, mBaseLineVariants, head);
      }

      if (head.matches()) {
        addIfBetter(head, sortedPaths);
      }
    }
    //System.err.println("Best: " + best);
    Diagnostic.userLog("Reference " + mTemplateName + " had maximum path complexity of " + maxPaths + " at " + maxPathsRegion);
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
    //System.err.println("Skip to next variant " + lastTemplatePos);
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
    final Map.Entry<Integer, Variant> entry = map.ceilingEntry(path.getPosition());
    return entry == null ? null : entry.getValue();
  }

  static <T extends Variant> TreeMap<Integer, Variant> buildMap(Collection<T> variations) {
    final TreeMap<Integer, Variant> map = new TreeMap<>();
    for (final Variant v : variations) {
      // TODO when you have a pure insert immediately prior to a snp/mnp, they end up with the same start position,
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

}
