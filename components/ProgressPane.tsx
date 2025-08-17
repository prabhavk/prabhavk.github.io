"use client";

import { useEffect, useRef } from "react";
import { useProgress } from "./ProgressProvider";

export default function ProgressPane() {
  const { jobId, status, logs, clear } = useProgress();
  const boxRef = useRef<HTMLDivElement>(null);

  // auto-scroll to bottom as logs arrive
  useEffect(() => {
    const el = boxRef.current;
    if (el) el.scrollTop = el.scrollHeight;
  }, [logs]);

  return (
    <aside className="w-full lg:w-80 xl:w-96 border-l bg-white p-4 sticky top-0 h-[calc(100vh-0px)]">
      <div className="flex items-center justify-between mb-2">
        <h2 className="font-semibold">Progress</h2>
        <span className="text-xs rounded px-2 py-0.5 bg-gray-100 text-gray-700">
          {status.toUpperCase()}
        </span>
      </div>
      <div className="text-xs text-gray-600 mb-2">Job: {jobId ?? "—"}</div>

      <div
        ref={boxRef}
        className="h-[70vh] overflow-auto rounded border bg-gray-50 p-2 font-mono text-xs text-black"
      >
        {logs.length === 0 ? <div className="text-gray-400">No output yet…</div> :
          logs.map((l, i) => <div key={i}>{l}</div>)}
      </div>

      <div className="mt-3 flex gap-2">
        <button
          onClick={() => navigator.clipboard.writeText(logs.join("\n"))}
          className="px-3 py-1 rounded bg-blue-600 text-white text-sm hover:bg-blue-700"
        >
          Copy
        </button>
        <button
          onClick={clear}
          className="px-3 py-1 rounded bg-gray-200 text-sm hover:bg-gray-300"
        >
          Clear
        </button>
      </div>
    </aside>
  );
}
