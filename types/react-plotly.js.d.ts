declare module "react-plotly.js" {
  import * as React from "react";

  export interface PlotParams {
    data: Partial<Plotly.Data>[];
    layout?: Partial<Plotly.Layout>;
    config?: Partial<Plotly.Config>;
    style?: React.CSSProperties;
    className?: string;
    onInitialized?: (figure: any, graphDiv: any) => void;
    onUpdate?: (figure: any, graphDiv: any) => void;
    onClick?: (event: any) => void;
    onHover?: (event: any) => void;
    onUnhover?: (event: any) => void;
  }

  export default class Plot extends React.Component<PlotParams> {}
}
