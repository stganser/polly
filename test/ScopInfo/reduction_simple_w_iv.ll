; RUN: opt %loadPolly -polly-scops -analyze < %s | FileCheck %s
;
; CHECK: Reduction like: 1
;
; void f(int* sum) {
;   for (int i = 0; i <= 100; i++)
;     sum += i * 3;
; }
target datalayout = "e-m:e-p:32:32-i64:64-v128:64:128-n32-S64"

define void @f(i32* %sum) {
entry:
  br label %for.cond

for.cond:                                         ; preds = %for.cond, %entr
  %i1.0 = phi i32 [ 0, %entry ], [ %inc, %for.cond ]
  %sum.reload = load i32* %sum
  %mul = mul nsw i32 %i1.0, 3
  %add = add nsw i32 %sum.reload, %mul
  %inc = add nsw i32 %i1.0, 1
  store i32 %add, i32* %sum
  %cmp = icmp slt i32 %i1.0, 100
  br i1 %cmp, label %for.cond, label %for.end

for.end:                                          ; preds = %for.cond
  ret void
}