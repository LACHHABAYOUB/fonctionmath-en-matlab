

	m0=1.00000597682;
	m1=9.54786104043e-4;
	m2=2.85583733151e-4;
	G=2.95912208286e-4;


	q0=[0;0;0; -3.5023653; -3.8169847; -1.5507963;9.0755314; -3.0458353; -1.6483708];
	p0=[0;0;0;0.00565429*m1;-0.00412490*m1;-0.00190589*m1;0.00168318*m2;0.00483525*m2;0.00192462*m2];



	[vp,vq]=euler(3000,50,p0,q0)
	


	vq(4:6,:)=vq(4:6,:)-vq(1:3,:);    
	vq(7:9,:)=vq(7:9,:)-vq(1:3,:);   
	plot([vq(4,:)’,vq(7,:)’],[vq(5,:)’,vq(8,:)’]);
	


	comet3d([vq(4,:)’,vq(7,:)’],[vq(5,:)’,vq(8,:)’],[vq(6,:)’,vq(9,:)’]);