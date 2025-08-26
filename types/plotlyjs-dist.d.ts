declare module "plotly.js-dist" {
  export type ImageFormat = "png" | "jpeg" | "svg" | "webp";

  export interface DownloadImageOptions {
    format?: ImageFormat;
    width?: number;
    height?: number;
    filename?: string;
  }

  export interface PlotlyLike {
    downloadImage: (gd: HTMLElement, opts: DownloadImageOptions) => Promise<string>;
  }

  const Plotly: PlotlyLike;
  export default Plotly;
}