import { Suspense } from "react";
import ResultsClient from "@/app/results/ResultsClient";

export default function ResultsPage() {
  return (
    <Suspense fallback={<div className="p-4">Loading results…</div>}>
      <ResultsClient />
    </Suspense>
  );
}
