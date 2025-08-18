"use client";

import Link from "next/link";
import { usePathname } from "next/navigation";

export default function Navbar() {
  const pathname = usePathname();

  const tab = (href: string, label: string) => {
    const active = pathname === href || pathname?.startsWith(href + "/");
    return (
      <Link
        href={href}
        className={`px-4 py-2 rounded-md text-sm font-semibold transition-colors ${
          active
            ? "bg-gray-200 text-black shadow-md"
            : "bg-gray-700 text-white hover:bg-yellow-500"
        }`}
      >
        {label}
      </Link>
    );
  };

  return (
    <nav className="w-full border-b bg-black">
      <div className="max-w-6xl mx-auto px-4 py-3 flex gap-3">
        {tab("/", "Input")}
        {tab("/wilcoxon", "Wilcoxon Rank Sum Tests & Violin Plots")}
        {tab("/flower", "Flower Plots")}
      </div>
    </nav>
  );
}
