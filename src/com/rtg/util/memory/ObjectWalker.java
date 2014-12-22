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
package com.rtg.util.memory;

import java.lang.ref.Reference;
import java.lang.reflect.Array;
import java.lang.reflect.Field;
import java.util.Iterator;
import java.util.Stack;

/**
 * Visits each of the objects reachable from a particular object.<p>
 *
 */
public class ObjectWalker {

  /** All objects encountered. Set<Object> */
  private IdentitySet mRefs = null;

  /** During traversal, contain path to the current object */
  private Stack<Object> mPath = null;

  /** Optional filter used to accept or reject objects to be walked */
  private ExcludeClassFilter mWalkFilter = null;

  /** Whether to follow Reference objects (e.g.: SoftReferences) */
  private boolean mFollowRefs = true;

  /**
   * Sets a filter that determines whether an object should be walked.
   * This determines whether the object will have its fields scanned
   * for walking.
   *
   * @param filter an <code>ExcludeClassFilter</code> value. May be null.
   */
  public void setWalkFilter(final ExcludeClassFilter filter) {
    mWalkFilter = filter;
  }


  /**
   * Sets whether Reference objects will be followed.  The default is
   * true.  For example, you can set this to false to not pass through
   * SoftReferences.
   *
   * @param follow a <code>boolean</code> value
   */
  public void setFollowReferences(final boolean follow) {
    mFollowRefs = follow;
  }


  /**
   * The method that gets called for each object visited.
   *
   * @param path a <code>Stack</code> containing the path of objects
   * walked to reach this object.
   * @param value the <code>Object</code> being visited.
   */
  protected void visitObject(final Stack<Object> path, final Object value) {
    if (path.isEmpty()) {
      System.err.println("Found object of type: " + value.getClass());
    } else {
      System.err.println("Found object of type: " + value.getClass() + " in " + path.peek().getClass());
    }
  }


  private void walkInternal(final Object obj) {
    if ((obj == null) || (mRefs.contains(obj))) {
      return;
    }
    mRefs.add(obj);
    visitObject(mPath, obj);
    if (!(mWalkFilter == null || mWalkFilter.accept(mPath, obj))) {
      return;
    }
    if (obj instanceof Reference) {
      if (!mFollowRefs) {
        return;
      } else {
        walkInternal(((Reference<?>) obj).get());
      }
    }
    mPath.push(obj);
    //get the class info and calculate size
    final Class<?> classId = obj.getClass();
    final ClassMemory.Info info = ClassMemory.getMemoryInfo(classId);
    if (classId.isArray()) {
      final int length = Array.getLength(obj);
      if (info.isArray()) {
        //but not primitive array
        for (int i = 0; i < length; i++) {
          walkInternal(Array.get(obj, i));
        }
      }
    } else {
      //scan all subfields containing Objects and arrays
      for (final Iterator<Field> iter = info.getNonprimitiveIterator(); iter.hasNext(); ) {
        final Field field = iter.next();
        try {
          walkInternal(field.get(obj));
        } catch (final IllegalAccessException e) {
          throw new RuntimeException("Illegally attempting to access value of field:"
                                     + field + " on object:" + obj + " of class:" + classId);
        }
      }
    }
    mPath.pop();
  }


  /**
   * Walks the object.
   *
   * @param obj The Object to walk.
   */
  public void walk(final Object obj) {
    if (obj != null) {
      mRefs = new IdentitySet();
      mPath = new Stack<>();
      walkInternal(obj);
      mPath = null;
      mRefs = null;
    }
  }

}
