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

import java.util.Arrays;
import java.util.Collection;
import java.util.TreeSet;

import com.reeltwo.jumble.annotations.TestClass;
import com.rtg.launcher.GlobalFlags;
import com.rtg.util.BasicLinkedListNode;
import com.rtg.util.diagnostic.Diagnostic;

/**
 * Find the path through the two sequences of variants that best reconciles them.
 */
@TestClass("com.rtg.variant.eval.PathTest")
public final class PathFinder {

  //private static final boolean DUMP = false;

  private static final int MAX_COMPLEXITY = GlobalFlags.getIntegerValue(GlobalFlags.VCFEVAL_MAX_PATHS); // Threshold on number of unresolved paths
  private static final int MAX_ITERATIONS = GlobalFlags.getIntegerValue(GlobalFlags.VCFEVAL_MAX_ITERATIONS);  // Threshold on number of iterations since last sync point

  private final byte[] mTemplate;
  private final String mTemplateName;
  private final Variant[] mCalledVariants;
  private final Variant[] mBaseLineVariants;

  private <T extends Variant> PathFinder(byte[] template, String templateName, Collection<T> calledVariants, Collection<T> baseLineVariants) {
    mTemplate = template;
    mTemplateName = templateName;

    mCalledVariants = calledVariants.toArray(new Variant[calledVariants.size()]);
    Arrays.sort(mCalledVariants, DetectedVariant.NATURAL_COMPARATOR);

    mBaseLineVariants = baseLineVariants.toArray(new Variant[baseLineVariants.size()]);
    Arrays.sort(mBaseLineVariants, DetectedVariant.NATURAL_COMPARATOR);
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

  Path bestPath() {
    // make it easy to find variants
    final TreeSet<Path> sortedPaths = new TreeSet<>();
    sortedPaths.add(new Path(mTemplate));
    Path best = null;
    int maxPaths = 0;
    String maxPathsRegion = "";
    int currentIterations = 0;
    int currentMaxIterations = 0;
    int currentMax = 0;
    int currentMaxPos = 0;
    Path lastSyncPath = null;
    int lastSyncPos = 0;
    String lastWarnMessage = null;
    while (sortedPaths.size() > 0) {
      currentMax = Math.max(currentMax, sortedPaths.size());
      currentMaxIterations = Math.max(currentMaxIterations, currentIterations++);
      Path head = sortedPaths.pollFirst();
      // if (DUMP) System.err.println("Size: " + (sortedPaths.size() + 1) + " Range:" + (lastSyncPos + 1) + "-" + (currentMaxPos + 1) + " LocalIterations: " + currentIterations + "\n\nHead: " + head);
      if (sortedPaths.size() == 0) { // Only one path currently in play
        if (lastWarnMessage != null) { // Issue a warning if we encountered problems during the previous region
          Diagnostic.warning(lastWarnMessage);
          lastWarnMessage = null;
        }
        final int currentSyncPos = head.mCalledPath.getPosition();
        if (currentMax > maxPaths) {
          maxPathsRegion = mTemplateName + ":" + (lastSyncPos + 1) + "-" + (currentSyncPos + 1);
          maxPaths = currentMax;
          Diagnostic.developerLog("Maximum path complexity now " + maxPaths + ", at " + maxPathsRegion + " with "  + currentIterations + " iterations");
        }
        currentMax = 0;
        currentIterations = 0;
        lastSyncPos = currentSyncPos;
        lastSyncPath = head;
      } else if (sortedPaths.size() > MAX_COMPLEXITY || currentIterations > MAX_ITERATIONS) {
        lastWarnMessage = "Evaluation too complex (" + sortedPaths.size() + " unresolved paths, " + currentIterations + " iterations) at reference region " + mTemplateName + ":" + (lastSyncPos + 1) + "-" + (currentMaxPos + 2) + ". Variants in this region will not be included in results.";
        sortedPaths.clear();    // Drop all paths currently in play
        currentIterations = 0;
        head = lastSyncPath;    // Create new head containing path up until last sync point
        skipVariantsTo(head.mCalledPath, mCalledVariants, currentMaxPos + 1);
        skipVariantsTo(head.mBaselinePath, mBaseLineVariants, currentMaxPos + 1);
      }
      if (head.finished()) {
        // Path is done. Remember the best
        final BasicLinkedListNode<Integer> syncPoints = new BasicLinkedListNode<>(head.mCalledPath.getPosition(), head.mSyncPointList);
        best = better(best, new Path(head, syncPoints));
        continue;
      }
      final int aVarIndex = nextVariant(head.mCalledPath, mCalledVariants);
      final Variant aVar = (aVarIndex == -1) ? null : mCalledVariants[aVarIndex];
      if (aVar != null && head.mCalledPath.isNew(aVar)) {
        currentMaxPos = Math.max(currentMaxPos, aVar.getStart());
        //Adding a new variant to A side
        // if (DUMP) System.err.println("Add alternatives to called " + aVar);
        addIfBetter(head.addAVariant(aVar, aVarIndex), sortedPaths);
        continue;
      }
      final int bVarIndex = nextVariant(head.mBaselinePath, mBaseLineVariants);
      final Variant bVar = (bVarIndex == -1) ? null : mBaseLineVariants[bVarIndex];
      if (bVar != null && head.mBaselinePath.isNew(bVar)) {
        currentMaxPos = Math.max(currentMaxPos, bVar.getStart());
        //Adding a new variant to B side
        // if (DUMP) System.err.println("Add alternatives to baseline " + bVar);
        addIfBetter(head.addBVariant(bVar, bVarIndex), sortedPaths);
        continue;
      }

      head.step();

      if (head.inSync()) {
        skipToNextVariant(head);
        // if (DUMP) System.err.println("In sync, skipping: " + head);
        // } else {
        // if (DUMP) System.err.println("Not in sync");
      }

      if (head.matches()) {
        // if (DUMP) System.err.println("Head matches, keeping");
        addIfBetter(head, sortedPaths);
        // } else {
        // if (DUMP) System.err.println("Head mismatch, discard");
      }
    }
    //System.err.println("Best: " + best);
    Diagnostic.userLog("Reference " + mTemplateName + " had maximum path complexity of " + maxPaths + " at " + maxPathsRegion);
    return best;
  }

  /**
   * Move the half path to the specified position, ignoring any intervening variants.
   * @param path the half path to skip
   * @param variants all variants
   * @param maxPos reference position to move to.
   */
  private void skipVariantsTo(HalfPath path, Variant[] variants, int maxPos) {
    int varIndex = path.getVariantIndex();
    while (varIndex < variants.length && variants[varIndex].getStart() < maxPos) {
      varIndex++;
    }
    varIndex--;
    Diagnostic.developerLog("Skipped to maxPos: " + maxPos + ". Variant index " + path.getVariantIndex() + " -> " + varIndex);
    path.setVariantIndex(varIndex);
    path.moveForward(maxPos);
  }


  /**
   * Move the path to just before the next variant
   * @param head the path to skip forward
   */
  private void skipToNextVariant(final Path head) {
    final int aNext = futureVariantPos(head.mCalledPath, mCalledVariants);
    final int bNext = futureVariantPos(head.mBaselinePath, mBaseLineVariants);
    final int lastTemplatePos = mTemplate.length - 1;
    //System.err.println("Skip to next variant " + lastTemplatePos);
    // -1 because we want to be before the position
    final int nextPos = Math.min(Math.min(aNext, bNext), lastTemplatePos) - 1;
    assert head.mCalledPath.getPosition() == head.mBaselinePath.getPosition();
    if (nextPos > head.mCalledPath.getPosition()) {
      head.moveForward(nextPos);
    }
  }

  /* Gets the next upstream variant position */
  int futureVariantPos(HalfPath path, Variant[] variants) {
    final int nextIdx = path.getVariantIndex() + 1;
    if (nextIdx >= variants.length) {
      return mTemplate.length - 1;
    } else {
      return variants[nextIdx].getStart();
    }
  }


  /* Gets the next variant that should be enqueued to the supplied HalfPath at the current position, or null if there is none  */
  static int nextVariant(HalfPath path, Variant[] variants) {
    final int nextIdx = path.getVariantIndex() + 1;
    if (nextIdx >= variants.length) {
      return -1;
    }
    final Variant nextVar = variants[nextIdx];

    final int startPos;
    if (path.wantsFutureVariantBases()) {
      startPos = Math.max(path.getVariantEndPosition(), path.getPosition() + 1);
    } else {
      startPos = path.getPosition() + 1;
    }
    return nextVar.getStart() == startPos ? nextIdx : -1;
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
