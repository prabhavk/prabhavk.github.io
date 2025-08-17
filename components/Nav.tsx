'use client';

import Link from "next/link";
import { usePathname } from "next/navigation";

export default function Nav() {
  const pathname = usePathname();
  const linkBase = "px-4 py-2 rounded-md text-sm font-medium";
  const active   = "bg-blue-600 text-white";
  const idle     = "bg-gray-100 text-black hover:bg-gray-200";

  return (
    <nav className="w-full flex gap-2 p-4 border-b bg-white justify-start">
      <Link href="/" className={`${linkBase} ${pathname === "/" ? active : idle}`}>
        Input
      </Link>
      <Link href="/wilcoxon" className={`${linkBase} ${pathname?.startsWith("/wilcoxon") ? active : idle}`}>
        Wilcoxon Rank Sum Tests & Violin Plots
      </Link>
      <Link href="/flower" className={`${linkBase} ${pathname?.startsWith("/flower") ? active : idle}`}>
        Flower Plots
      </Link>
    </nav>
  );
}
