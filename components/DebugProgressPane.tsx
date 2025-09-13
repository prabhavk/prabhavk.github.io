// components/DebugProgressPane.tsx
"use client";

import React from "react";
import { useShallow } from "zustand/react/shallow";
import { useProgressStore } from "@/stores/progress";

export default function DebugProgressPane({
  height = 500,
  className = "",
}: {
  height?: number;
  className?: string;
}) {
  const [
    logs,
    clearLogs,
    currentMethod,
    currentRoot,
    currentRep,
    totalReps,
    nodeIndex,
    nodeTotal,
    lastLayer,
  ] = useProgressStore(
    useShallow((s) => [
      s.logs,            // NEW
      s.clearLogs,       // NEW
      s.currentMethod,
      s.currentRoot,
      s.currentRep,
      s.totalReps,
      s.nodeIndex,
      s.nodeTotal,
      s.lastLayer,
    ] as const)
  );

  const endRef = React.useRef<HTMLDivElement>(null);
  React.useEffect(() => {
    endRef.current?.scrollIntoView({ behavior: "smooth" });
  }, [logs.length]);

  return (
    <div
      className={`rounded-2xl border border-stone-300 bg-white text-black p-3 ${className}`}
      style={{ height }}
    >
      <div className="flex items-center justify-between mb-2">
        <div className="text-sm">
          <div className="font-semibold">Progress</div>
          <div className="text-xs text-black/70">
            rep {currentRep ?? 0}/{totalReps ?? 0} · method {currentMethod ?? "—"} · layer{" "}
            {lastLayer ?? "—"} · node {nodeIndex ?? 0}/{nodeTotal ?? 0}
            {currentRoot ? ` (${currentRoot})` : ""}
          </div>
        </div>

        <button
          onClick={clearLogs}
          className="text-xs border px-2 py-1 rounded hover:bg-stone-50"
          disabled={logs.length === 0}
        >
          Clear
        </button>
      </div>

      <div className="h-[calc(100%-2.5rem)] overflow-auto rounded bg-stone-50 border border-stone-200 p-2">
        <pre className="text-xs leading-5 whitespace-pre-wrap break-words">
          {logs.length === 0 ? (
            <span className="text-black/50">
              Waiting for logs… (enable “Upsert UI logs” or ensure the worker is running)
            </span>
          ) : (
            logs.map((l, i) => <div key={i}>{l}</div>)
          )}
          <div ref={endRef} />
        </pre>
      </div>
    </div>
  );
}
