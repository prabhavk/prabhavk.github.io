// hooks/useEmtrProgress.ts
import { useEffect, useMemo, useRef, useState } from "react";

export type Method = "parsimony" | "dirichlet" | "hss";

type EmtrState = {
  currentMethod: Method | null;
  currentNode: string | null;
  
  completedReps: number;
  totalReps: number;
  
  nodeOrder: Record<Method, string[]>;
};

type Options = {
  worker: Worker | null;
  totalReps: number;
};

export function useEmtrProgress({ worker, totalReps }: Options): EmtrState {
  const [currentMethod, setCurrentMethod] = useState<Method | null>(null);
  const [currentNode, setCurrentNode] = useState<string | null>(null);
  const [completedReps, setCompletedReps] = useState(0);
  const nodeOrderRef = useRef<Record<Method, string[]>>({
    parsimony: [],
    dirichlet: [],
    hss: [],
  });
  const seenRepForMethod = useRef<Record<Method, Set<number>>>({
    parsimony: new Set(),
    dirichlet: new Set(),
    hss: new Set(),
  });

  useEffect(() => {
    if (!worker) return;
    const onMsg = (evt: MessageEvent) => {
      const { data } = evt;
      if (!data) return;

      
      if (data.type === "emtrifle" && data.payload) {
        const p = data.payload;
        const method: Method = (String(p.method || "") as Method) || "parsimony";
        const root = String(p.root || "");
        const rep = Number(p.rep) || 0;

        setCurrentMethod(method);
        if (root) setCurrentNode(root);

      
        const arr = nodeOrderRef.current[method];
        if (root && !arr.includes(root)) arr.push(root);

      
        if (rep > 0) {
          const bucket = seenRepForMethod.current[method];
          const before = bucket.size;
          bucket.add(rep);
      
          const totalSeen = new Set<number>([
            ...seenRepForMethod.current.parsimony,
            ...seenRepForMethod.current.dirichlet,
            ...seenRepForMethod.current.hss,
          ]).size;
          if (totalSeen !== completedReps) setCompletedReps(totalSeen);
        }
      }
    };
    worker.addEventListener("message", onMsg);
    return () => worker.removeEventListener("message", onMsg);
  }, [worker, completedReps]);

  const nodeOrder = useMemo(() => ({ ...nodeOrderRef.current }), []);
  return { currentMethod, currentNode, completedReps, totalReps, nodeOrder };
}
