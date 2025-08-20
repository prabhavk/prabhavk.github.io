// app/(input)/layout.tsx
import ProgressProvider from "@/components/ProgressProvider";
import ProgressPane from "@/components/ProgressPane";

export default function InputLayout({ children }: { children: React.ReactNode }) {
  return (
    <ProgressProvider>
      {/* Full-bleed wrapper so the pane hugs the right edge */}
      <div className="relative left-1/2 right-1/2 -ml-[50vw] -mr-[50vw] w-screen">
        {/* Add left gutter; keep right flush */}
        <div className="grid grid-cols-1 lg:grid-cols-[minmax(0,3fr)_minmax(0,2fr)] gap-3
                        pr-0 pl-6 sm:pl-8 lg:pl-10 xl:pl-12 2xl:pl-16">
          <main id="main" className="min-w-0">{children}</main>
          <aside className="lg:sticky lg:top-4 self-start max-h-[calc(100vh-1rem)] overflow-auto">
            <ProgressPane />
          </aside>
        </div>
      </div>
    </ProgressProvider>
  );
}
