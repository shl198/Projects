import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import java.io.File
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode.GVCF
import org.broadinstitute.gatk.utils.variant.GATKVCFIndexType.LINEAR


class HaplotypeCallerDNAgVCF extends QScript {
	
	// Define the input
	
	@Input(doc="The reference genome file organism.fa", shortName="R", required=true)
	var refgenomeFile: File = _

	@Input(doc="Input recalibration file name.recal.bam", shortName="I", required=true)
	var recaliFile: File = _

	@Output(doc="Output gvcf file name.gvcf", shortName="O", required=true)
	var gvcfFile: File = _

	@Argument(doc="index paramter", shortName="variant_index_parameter", required=false)
	var index_parameter: Int= 128000

	@Argument(doc="number of scatter", shortName="sg", required=false)
	var sg: Int= _
	

	def script() {
		val haplotypeCaller = new HaplotypeCaller
		haplotypeCaller.reference_sequence = new File(refgenomeFile)
		haplotypeCaller.input_file :+= new File(recaliFile)
		haplotypeCaller.emitRefConfidence = GVCF
		haplotypeCaller.variant_index_type = LINEAR
		haplotypeCaller.variant_index_parameter = index_parameter
		haplotypeCaller.out = new File(gvcfFile)
		haplotypeCaller.scatterCount = sg
		add(haplotypeCaller)

	}
}
